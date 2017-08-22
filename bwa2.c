#include <inttypes.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "hash.h"
#include "bwtaln.h"
#include "bwt.h"
#include "aryana_args.h"
#include "bwa2.h"
#include "aligner.h"
#include "bntseq.h"
#include "utils.h"
#include "const.h"

//#include "aligner.h"
#include "sam.h"
#define MAX(a, b) (a > b) ? a : b
char atom2[4]= {'A','C','G','T'};
int counter = 0;
double base=10.0;
double scmax=100.0;
bwa_seq_t *seqs;
int n_seqs = 0;
int tot_seqs = 0;
static int output_buffer_warning = 0;
pthread_mutex_t input = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t output = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int tid;
    const gap_opt_t *opt;
    bwt_t * bwt;
    aryana_args *args;
    bwa_seqio_t * ks, *ks2;
    uint64_t * reference, reference_size, reference_reminder;
    bwtint_t * offset, offInd;
} global_vars;

// The Data structure used for preserving the initial order of the reads
typedef struct {
    char * str;
    unsigned long long read_num;
} outBuf_t;
outBuf_t * outBuf;
int running_threads_num;
unsigned long long outBufIndex, outBufLen;

void compute_penalty(aryana_args * args, penalty_t * p) {
    if (args->mismatch_limit >= 0 && p->mismatch_num > args->mismatch_limit) p->penalty = maxPenalty;
    else p->penalty = args->mismatch_penalty * p->mismatch_num + args->gap_open_penalty * p->gap_open_num + args->gap_ext_penalty * p->gap_ext_num;
}

void OutputBufferError() {
    fprintf(stderr, "Error: Output buffer size is insufficient for preserving the order of reads.\n");
    fprintf(stderr, "Either set a larger size for output buffer using -B/--buffer argument, or don't use -O/--order argument.\n");
    exit(1);
}

void OutputBufferWarning() {
    if (output_buffer_warning) return;
    output_buffer_warning = 1;
    fprintf(stderr, "Warning: Output buffer size is low for preserving the order of reads. This might reduce the speed.\n");
    fprintf(stderr, "Set a larger size for output buffer using -B/--buffer argument to omit this warning.\n");
}

void PrintSeq(unsigned char * seq, int len, int endl) {
    int i;
    for (i = 0; i < len; i++) fprintf(stderr, "%c", atom2[seq[i]]);
    if (endl) fprintf(stderr, "\n");
}

void PrintRefSeq(uint64_t * reference, unsigned long long start, unsigned long long end, unsigned long long seq_len, int endl) {
    unsigned long long k;
    if (end >= seq_len) end = seq_len - 1;
    for (k = start; k <= end; k++)
        fprintf(stderr, "%c", atom2[(unsigned short)getNuc(k, reference, seq_len)]);
    if (endl) fprintf(stderr, "\n");
}

void GetRefSeq(uint64_t * reference, unsigned long long start, unsigned long long end, unsigned long long seq_len, char * seq) {
    unsigned long long k;
    for (k = start; k <= end; k++)
        seq[k-start] = atom2[(unsigned short)getNuc(k, reference, seq_len)];
    seq[k-start] = 0;
}

void compute_cigar_penalties(global_vars *g, int candidates_num, int candidates_size, int * candidates, char * cigar[], bwa_seq_t * seq, hash_element * table, int **d, char **arr, char *tmp_cigar, penalty_t * penalty) {
    int i;
    for (i= candidates_num; i< candidates_size; i++) // The array is started to fill from the last.
        if (candidates[i]!=-1)
        {
            create_cigar(g->args, &table[candidates[i]], cigar[i], seq->len, seq->seq, seq->qual, g->bwt->seq_len,d,arr, tmp_cigar, penalty + i, g->reference, g->args->ignore);
            compute_penalty(g->args, penalty + i);
            seq->index=table[candidates[i]].index;
            if (g->args->debug >= 1) {
                fprintf(stderr, "index: %llu, cigar: %s, mimatch %d, gap open %d, gap ext %d, penalty %d\n", (unsigned long long) table[candidates[i]].index,  cigar[i], penalty[i].mismatch_num, penalty[i].gap_open_num, penalty[i].gap_ext_num, penalty[i].penalty);
                PrintSeq(seq->seq, seq->len, 1);
                PrintRefSeq(g->reference, table[candidates[i]].index, table[candidates[i]].index + seq->len - 1, g->bwt->seq_len, 1);
            }
        }
}

int valid_pair(global_vars * g, hash_element *hit1, hash_element * hit2, bwtint_t len1, bwtint_t len2) {
    bwtint_t seqlen = g->bwt->seq_len;
    bwtint_t index[2], len[2], fragLen;
    index[0] = hit1->index;
    index[1] = hit2->index;
    len[0] = len1;
    len[1] = len2;
    int rev[2], i;
    bzero(rev, sizeof(rev));
    for (i = 0; i < 2; i++)
        if (index[i] > seqlen / 2) {
            index[i] = (seqlen / 2) - (index[i] - (seqlen / 2)) - len[i];
            rev[i] = 1;
        }
    if (index[1] > index[0]) fragLen = index[1] - index[0] + len[1];
    else fragLen = index[0] - index[1] + len[0];
    if (fragLen < g->args->min_dis || fragLen > g->args->max_dis) return 0; 	// The length of fragment does not match the input arguments
    if (get_chr(index[0], seqlen, g->offset, g->offInd) != get_chr(index[1], seqlen, g->offset, g->offInd)) return 0;
    switch (g->args->orientation) {
    case orien_all:
        break;
    case orien_ff:
        if (rev[0] != rev[1]) return 0;
        break;
    case orien_fr:
        if (rev[0] == rev[1] || (rev[0] && index[0] < index[1]) || (rev[1] && index[1] < index[0])) return 0;
        break;
    case orien_rf:
        if (rev[0] == rev[1] || (rev[1] && index[0] < index[1]) || (rev[0] && index[1] < index[0])) return 0;
        break;
    }
    return 1;	// Valid pair
}

void compute_mapq(global_vars * g, int * can, int * can2, int num, penalty_t * p, penalty_t * p2, hash_element * t, hash_element * t2, char * c[], char * c2[], int paired) {
    /*int i;
    for (i = 0; i < num; i++) {
        p[i].mapq = p[i].mapq/;
        if (paired) p2[i].mapq = MAX(1, 41 - p2[i].penalty/5);
    }*/
}

void find_best_candidates(global_vars * g, int candidates_num, int candidates_num2, int candidates_size, int * candidates, int * candidates2, penalty_t * penalty, penalty_t * penalty2, int * best_candidates, int * best_candidates2, int * best_num, hash_element *table, hash_element * table2, bwtint_t len, bwtint_t len2,  char * cigar[], char * cigar2[], int paired) {
    *best_num = 0;
    int min_penalty = maxPenalty, i, j;
    for (i= candidates_num; i< candidates_size; i++) // The array is started to fill from the last.
        if (candidates[i]!=-1) {
            if (! paired) {	// Single-end reads
                if (min_penalty > penalty[i].penalty) // improving single-end alignment penalty
                {
                    min_penalty = penalty[i].penalty;
                    *best_num = 0;
                }
                if (min_penalty < maxPenalty && min_penalty  == penalty[i].penalty) {
                    best_candidates[*best_num] = candidates[i]; // a new single-end alignment with the same penalty
                    penalty[*best_num] = penalty[i];
                    strcpy(cigar[*best_num], cigar[i]);
                    (*best_num)++;
                }
            } else { // Paired-end reads
                char tcigar[g->args->potents][MAX_CIGAR_SIZE], tcigar2[g->args->potents][MAX_CIGAR_SIZE];
                penalty_t tpenalty[g->args->potents], tpenalty2[g->args->potents];
                for (j = candidates_num2; j < candidates_size; j++)
                {
                    if (candidates2[j] != -1 && valid_pair(g, table + candidates[i], table2 + candidates2[j], len, len2)) {
                        if (min_penalty > penalty[i].penalty + penalty2[j].penalty) { // improving paired-alignment penalty
                            min_penalty =  penalty[i].penalty + penalty2[j].penalty;
                            *best_num = 0;
                        }
                        if (*best_num < candidates_size && min_penalty < maxPenalty && min_penalty ==  penalty[i].penalty + penalty2[j].penalty) { // A new paired alignment with the same penalty
                            best_candidates[*best_num] = candidates[i];
                            best_candidates2[*best_num] = candidates2[j];
                            tpenalty[*best_num] = penalty[i];
                            tpenalty2[*best_num] = penalty2[j];
                            strcpy(tcigar[*best_num], cigar[i]);
                            strcpy(tcigar2[*best_num], cigar2[j]);
                            (*best_num)++;
                        }
                    }
                }
                for (j = 0; j < *best_num; j++) {
                    penalty[j] = tpenalty[j];
                    penalty2[j] = tpenalty2[j];
                    strcpy(cigar[j], tcigar[j]);
                    strcpy(cigar2[j], tcigar2[j]);
                }
            }
        }
    if (*best_num > 0)
        compute_mapq(g, best_candidates, best_candidates2, *best_num, penalty, penalty2, table, table2, cigar, cigar2, paired);
}

long long report_alignment_results(global_vars *g, char * buffer, char *cigar[], char * cigar2[], bwa_seq_t *seq, bwa_seq_t *seq2, hash_element *table, hash_element * table2, int * can, int * can2, int canNum, int canNum2, penalty_t * penalty, penalty_t * penalty2) {
    long long result = 0;
    int i, p, flag;
    if (!g->args->paired) {
        if (! canNum) return sam_generator(buffer + result, seq->name, sf_unmapped, 0, 0, 0, 0, 0, seq->seq, seq->qual, seq->len, 0, g->bwt, g->offset, g->offInd);
        p = rand() % canNum; // primary alignment
        result += sam_generator(buffer + result, seq->name, 0, penalty[p].mapq, table[can[p]].index, cigar[p], 0, 0, seq->seq, seq->qual, seq->len, 0, g->bwt, g->offset, g->offInd);
        if (g->args->report_multi)
            for (i = 0; i < canNum; i++)
                if (i != p)
                    result += sam_generator(buffer + result, seq->name, sf_not_primary, penalty[i].mapq, table[can[i]].index, cigar[i], 0, 0, seq->seq, seq->qual, seq->len, 0, g->bwt, g->offset, g->offInd);
        return result;
    }
    // Reporting the results for paired alignments
    flag = sf_paired;
    if (! canNum) {
        flag |= sf_unmapped + sf_mate_unmapped;
        result += sam_generator(buffer + result, seq->name, flag | sf_first_mate, 0, 0, 0, 0, 0, seq->seq, seq->qual, seq->len, 0, g->bwt, g->offset, g->offInd);
        result += sam_generator(buffer + result, seq2->name, flag | sf_second_mate, 0, 0, 0, 0, 0, seq2->seq, seq2->qual, seq2->len, 0, g->bwt, g->offset, g->offInd);
        return result;
    }
    flag |= sf_pair_mapped;
    p = rand() % canNum;
    result += sam_generator(buffer + result, seq->name, flag | sf_first_mate, penalty[p].mapq, table[can[p]].index, cigar[p], table2[can2[p]].index, cigar2[p], seq->seq, seq->qual, seq->len, seq2->len, g->bwt, g->offset, g->offInd);
    result += sam_generator(buffer + result, seq2->name, flag | sf_second_mate, penalty2[p].mapq, table2[can2[p]].index, cigar2[p], table[can[p]].index, cigar[p], seq2->seq, seq2->qual, seq2->len, seq->len, g->bwt, g->offset, g->offInd);
    if (g->args->report_multi) {
        flag |= sf_not_primary;
        for (i = 0; i < canNum; i++)
            if (i != p) {
                result += sam_generator(buffer + result, seq->name, flag | sf_first_mate, penalty[i].mapq, table[can[i]].index, cigar[i], table2[can2[i]].index, cigar2[i], seq->seq, seq->qual, seq->len, seq2->len, g->bwt, g->offset, g->offInd);
                result += sam_generator(buffer + result, seq2->name, flag | sf_second_mate, penalty2[i].mapq, table2[can2[i]].index, cigar2[i], table[can[i]].index, cigar[i], seq2->seq, seq2->qual, seq2->len, seq->len, g->bwt, g->offset, g->offInd);
            }
    }
    return result;
}

long long report_discordant_alignment_results(global_vars *g, char * buffer, char *cigar[], char * cigar2[], bwa_seq_t *seq, bwa_seq_t *seq2, hash_element *table, hash_element * table2, int * can, int * can2, int canNum, int canNum2, penalty_t * penalty, penalty_t * penalty2) {
    long long result = 0;
    int i, p, p2, flag;
    bwtint_t mate_ind;
    char * mate_cigar;
    if (canNum) p = rand() % canNum;
    if (canNum2) p2 = rand() % canNum2;
// Reporting matches for the first read
    flag = sf_paired | sf_first_mate;
    if (!canNum2) {
        flag |= sf_mate_unmapped;
        mate_ind = 0;
        mate_cigar = 0;
    } else {
        mate_ind = table2[can2[p2]].index;
        mate_cigar = cigar2[p2];
    }
    if (! canNum)
        result += sam_generator(buffer + result, seq->name, flag | sf_unmapped, 0, 0, 0, mate_ind, mate_cigar, seq->seq, seq->qual, seq->len, seq2->len, g->bwt, g->offset, g->offInd);
    else {
        result += sam_generator(buffer + result, seq->name, flag, penalty[p].mapq, table[can[p]].index, cigar[p], mate_ind, mate_cigar, seq->seq, seq->qual, seq->len, seq2->len, g->bwt, g->offset, g->offInd);
        if (g->args->report_multi) {
            flag |= sf_not_primary;
            for (i = 0; i < canNum; i++)
                if (i != p)
                    result += sam_generator(buffer + result, seq->name, flag, penalty[i].mapq, table[can[i]].index, cigar[i], mate_ind, mate_cigar, seq->seq, seq->qual, seq->len, seq2->len, g->bwt, g->offset, g->offInd);
        }
    }
// Reporting matches for the second read
    flag = sf_paired | sf_second_mate;
    if (!canNum) {
        flag |= sf_mate_unmapped;
        mate_ind = 0;
        mate_cigar = 0;
    } else {
        mate_ind = table[can[p]].index;
        mate_cigar = cigar[p];
    }
    if (! canNum2)
        result += sam_generator(buffer + result, seq2->name, flag | sf_unmapped, 0, 0, 0, mate_ind, mate_cigar, seq2->seq, seq2->qual, seq2->len, seq->len, g->bwt, g->offset, g->offInd);
    else {
        result += sam_generator(buffer + result, seq2->name, flag, penalty2[p2].mapq, table2[can2[p2]].index, cigar2[p2], mate_ind, mate_cigar, seq2->seq, seq2->qual, seq2->len, seq->len, g->bwt, g->offset, g->offInd);
        if (g->args->report_multi) {
            flag |= sf_not_primary;
            for (i = 0; i < canNum2; i++)
                if (i != p2)
                    result += sam_generator(buffer + result, seq2->name, flag, penalty2[i].mapq, table2[can2[i]].index, cigar2[i], mate_ind, mate_cigar, seq2->seq, seq2->qual, seq2->len, seq->len, g->bwt, g->offset, g->offInd);
        }
    }
    return result;
}

// This function aligns one (single or paired) read, finds the best alignment positions, computes the CIGAR sequence, and one or multiple lines of the SAM file for reporting the read alignments
// arr and d are the arrays used for smith_waterman dynamic programming. The reason to pass their pointers is to avoid creating them for every read
long long align_read(global_vars * g, char * buffer, char *cigar[], char * cigar2[], bwa_seq_t *seq, bwa_seq_t *seq2, hash_element *table, hash_element * table2, uint64_t level, int **d, char **arr, char* tmp_cigar) {
    buffer[0]='\0';
    int candidates_size=g->args->potents;
    int candidates[candidates_size], candidates2[candidates_size], candidates_num = candidates_size, candidates_num2 = candidates_size;
    int best_candidates[candidates_size], best_candidates2[candidates_size], best_num = 0, best_num2 = 0;
    penalty_t penalty[candidates_size], penalty2[candidates_size];
    int i=0;
    for (i=0; i< candidates_size; i++)
    {
        candidates[i]=-1;
        cigar[i][0]=0;
        if (g->args->paired) {
            candidates2[i] = -1;
            cigar2[i][0] = 0;
        }
    }
    //aligner
    aligner(g->bwt, seq->len, seq->seq, level, table, candidates, candidates_size, &candidates_num, g->args, g->reference);
    compute_cigar_penalties(g, candidates_num, candidates_size, candidates, cigar, seq, table, d, arr, tmp_cigar, penalty);
    if (g->args->paired) {
        aligner(g->bwt, seq2->len, seq2->seq,  level, table2, candidates2, candidates_size, &candidates_num2, g->args, g->reference);
        compute_cigar_penalties(g, candidates_num2, candidates_size, candidates2, cigar2, seq2, table2, d, arr, tmp_cigar, penalty2);
    }

    find_best_candidates(g, candidates_num, candidates_num2, candidates_size, candidates, candidates2, penalty, penalty2, best_candidates, best_candidates2, &best_num, table, table2, seq->len,  (g->args->paired) ? seq2->len : 0, cigar, cigar2, g->args->paired);
    if (g->args->paired && g->args->discordant && best_num == 0) { // No valid paired alignments but we should report "discordant" alignments: the best alignment for each read of the pair separately
        find_best_candidates(g, candidates_num, 0, candidates_size, candidates, 0, penalty, 0, best_candidates, 0, &best_num, table, 0, seq->len, 0, cigar, 0, 0);
        find_best_candidates(g, candidates_num2, 0, candidates_size, candidates2, 0, penalty2, 0, best_candidates2, 0, &best_num2, table2, 0, seq2->len, 0, cigar2, 0, 0);
        return report_discordant_alignment_results(g, buffer, cigar, cigar2, seq, seq2, table, table2, best_candidates, best_candidates2, best_num, best_num2, penalty, penalty2);
    }
    return report_alignment_results(g, buffer, cigar, cigar2, seq, seq2, table, table2, best_candidates, best_candidates2, best_num, best_num2, penalty, penalty2);
}



// This is the central function of each thread. It reads a group of reads, aligns each one, and puts the output to be printed in SAM file

void multiAligner(global_vars * g) {

    // Allocating memory and initializing variables
    bwa_seq_t *seqs, *seqs2 = 0;
    int n_seqs = 0, n_seqs2;
    int total_seqs = 0;
    bwtint_t j = 0;
    hash_element *table, *table2 = 0;

    table=(hash_element *)malloc(TABLESIZE*(sizeof (hash_element)));
    reset_hash(table);

    int **d=(int **)malloc(MAX_READ_SIZE*(sizeof (int *)));
    char **arr=(char **)malloc(MAX_READ_SIZE*(sizeof (char *)));
    char *tmp_cigar=(char *)malloc(MAX_READ_SIZE*(sizeof (char)));
    for (j=0; j<MAX_READ_SIZE; j++) {
        d[j]=(int *)malloc(300*(sizeof (int)));
        arr[j]=(char *)malloc(300*(sizeof (char)));
    }
    //seqs = (bwa_seq_t*)calloc(seq_num_per_read, sizeof(bwa_seq_t));
    char buffer[100*max_sam_line*(sizeof (char))];
    buffer[0] = '\0';
    long long lasttmpsize = 0;
    char **cigar = (char **) malloc(g->args->potents * sizeof(char *)), **cigar2 = 0;
    for (j=0; j < g->args->potents; j++)
        cigar[j]=(char *) malloc(MAX_CIGAR_SIZE*(sizeof (char)));

    if (g->args->paired) {
        table2=(hash_element *)malloc(TABLESIZE*(sizeof (hash_element)));
        reset_hash(table2);
        //seqs2 = (bwa_seq_t*) malloc(seq_num_per_read * sizeof(bwa_seq_t));
        //fprintf(stderr, "1, %d\n", seqs2);
        cigar2 =  (char **) malloc(g->args->potents * sizeof(char *));
        for (j=0; j < g->args->potents; j++)
            cigar2[j]=(char *) malloc(MAX_CIGAR_SIZE*(sizeof (char)));
    }

//    fprintf(stderr, "Thread %d starting...\n", g->tid);

    // The thread body
    while(true) {
        pthread_mutex_lock(&input);
        int read = 1;
        if((seqs = bwa_read_seq(g->ks, seq_num_per_read, &n_seqs, g->opt->mode, g->opt->trim_qual)) == 0 ||
                (g->args->paired && (seqs2 = bwa_read_seq(g->ks2, seq_num_per_read, &n_seqs2, g->opt->mode, g->opt->trim_qual)) == 0)) read = 0;
        pthread_mutex_unlock(&input);
        if (! read) break;
        if (g->args->paired && n_seqs != n_seqs2) {
            fprintf(stderr, "Error: The number of reads from fastq1 is not equal to fastq2\n");
            return;
        }
        lasttmpsize = 0;
        for(j=0; j<n_seqs; j++) {
            lasttmpsize+=align_read(g, buffer+lasttmpsize, cigar, cigar2, &seqs[j], (g->args->paired)? &seqs2[j] : 0, table, table2, total_seqs + j + 1, d,arr, tmp_cigar);

            if (g->args->order) { // Preserving order of reads
                unsigned long long read_num = seqs[j].read_num;
                int pos = read_num % outBufLen;
                while (outBuf[pos].read_num != 0) {
                    OutputBufferWarning();
                    usleep(1);
                }
                outBuf[pos].str = strdup(buffer);
                outBuf[pos].read_num = read_num; // Should be strictly after previous command of outBuf[pos].str assignment, for synchronization purposes.
                lasttmpsize = 0;
            } else {
                if(j % 100 == 0 || j == n_seqs - 1) {
                    pthread_mutex_lock(&output);
                    fputs(buffer, stdout);
                    pthread_mutex_unlock(&output);
                    lasttmpsize = 0;
                }
            }
        }
        total_seqs += n_seqs;
        bwa_free_read_seq(n_seqs, seqs);
        if (g->args->paired) bwa_free_read_seq(n_seqs, seqs2);
    }
    if (lasttmpsize > 0) {
        pthread_mutex_lock(&output);
        fputs(buffer, stdout);
        pthread_mutex_unlock(&output);
        lasttmpsize = 0;
    }
    // Free variables
    //free(seqs);
    for (j=0; j<MAX_READ_SIZE; j++) {
        free(d[j]);
        free(arr[j]);
    }
    for (j = 0; j < g->args->potents; j++)
        free(cigar[j]);
    free(cigar);
    free(table);
    if (g->args->paired) {
        free(table2);
        //fprintf(stderr, "2, %d\n", seqs2);
        //free(seqs2);
        for (j = 0; j < g->args->potents; j++)
            free(cigar2[j]);
        free(cigar2);
    }
    free(d);
    free(arr);
    free(tmp_cigar);
//    fprintf(stderr, "Thread %d finished.\n", g->tid);
}

void *worker2(void *data) {
    global_vars *g = (global_vars*)data;
    multiAligner(g);
    pthread_mutex_lock(&input);
    running_threads_num--;
    pthread_mutex_unlock(&input);
    return 0;
}

int ref_read(char * file_name, uint64_t ** reference_p, uint64_t * reference_size_p, uint64_t * reference_reminder_p) {
    struct stat file_info;
    if(stat(file_name , &file_info) == -1) {
        fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", file_name);
        return -1;
    }
    int fd = open(file_name , O_RDONLY);
    //FILE *fd;
    //fd = xopen(file_name, "rb");
    if(fd == -1) {
        fprintf(stderr, "Could not open the file %s\nplease make sure the file exists\n", file_name);
        return -1;
    }
    off_t file_size_bytes = file_info.st_size;
    *reference_size_p = ceil ( ((double)file_size_bytes) / (double)(sizeof(uint64_t)) );
    *reference_reminder_p = file_size_bytes % sizeof(uint64_t) ;
    *reference_p = (uint64_t *)malloc((*reference_size_p) * sizeof(uint64_t));
    //memset ( reference , 0 , reference_size * sizeof(uint64_t) ); Is this required?
    size_t read_size2 = 0;//there is a read_size defined above
    size_t signal;
    size_t total_size = (file_size_bytes);
    unsigned char *uc_buffer = (unsigned char *)(* reference_p);
    int counter=0;

    do {
        signal = read ( fd , (void *)uc_buffer , total_size - read_size2 );
        //signal = fread((void *)uc_buffer, )
        if ( signal == -1 )
        {
            fprintf(stderr, "Error: while writing to file\n");
            if ( close(fd) == -1 )
                fprintf(stderr, "Error: while closing file\n");
            return -1;
        }
        counter++;
        read_size2 += signal;
        uc_buffer += signal;
    }
    while ( read_size2 < total_size );
    if ( close(fd) == -1 )
    {
        fprintf(stderr, "Unable to close the file\n");
        return -1;
    }

    return 0;
}

char getNuc(uint64_t place, uint64_t * reference, uint64_t seq_len) {
    int rev=0;
    if(place >= (seq_len / 2))
    {
        place = (seq_len / 2) - (place - (seq_len / 2))-1;
        rev=1;
    }
    uint64_t block=place/(sizeof(bwtint_t)*4);
    int offset=place%(sizeof(bwtint_t)*4);
    uint64_t mask=3;
    mask=mask & (reference[block] >> (2*offset));
    if (rev==1)
        mask=3-mask;
    return mask;//atom[mask];
}

//void bwa_aln_core2(const char *prefix, const char *fn_fa, const char *fn_fa1, const char *fn_fa2, const //gap_opt_t *opt, pair_opt *args)
void bwa_aln_core2(aryana_args *args)
{
    bwt_t *bwt;
    bwa_seqio_t *ks, *ks2 = 0;
    gap_opt_t opt_tmp; //Afsoon: What's this?
    opt_tmp.mode = 0;
    opt_tmp.trim_qual = 0;
    gap_opt_t *opt = &opt_tmp;
    char *prefix = args->reference;
    char *fn_fa = args->read_file;
    char *fn_fa2 = args->read_file2;
    int i, j;
    uint64_t * reference;
    uint64_t reference_size, reference_reminder;
    bwtint_t *offset = (bwtint_t *)calloc(max_chrom_num, sizeof(bwtint_t));
    bwtint_t offInd = 0;
    ks = bwa_open_reads(opt->mode, fn_fa);
    if (args->paired)
        ks2 = bwa_open_reads(opt->mode, fn_fa2);

    // initialization

    // load BWT
    char *str = (char*)calloc(strlen(prefix) + 10, 1);
    strcpy(str, prefix);
    strcat(str, ".bwt");
    bwt = bwt_restore_bwt(str);
    strcpy(str, prefix);
    strcat(str, ".sa");
    bwt_restore_sa(str, bwt);
    strcpy(str, prefix);
    strcat(str, ".bin");
    ref_read(str, &reference, &reference_size, &reference_reminder);

    memset(offset, 0, max_chrom_num * sizeof(bwtint_t));
    strcpy(str, prefix);
    strcat(str, ".ann");
    FILE * ann = fopen(str, "r");
    free(str);

    char line[1000];
    if(fgets(line, sizeof line, ann) == NULL) {
        fprintf(stderr, "Error: Empty file\n");
        bwt_destroy(bwt);
        bwa_seq_close(ks);
        if (args->paired)
            bwa_seq_close(ks2);
        return;
    }
    fscanf(ann, "%d %s", &i, name[0]);
    while(fgets(line, sizeof line, ann) != NULL) {
        fscanf(ann, "%llu", (unsigned long long *) &offset[offInd++]);
        if(fgets(line, sizeof line, ann) == NULL) {
            fprintf(stderr, "Error: Empty file\n");
            bwt_destroy(bwt);
            bwa_seq_close(ks);
            if (args->paired)
                bwa_seq_close(ks2);
            return;
        }
        fscanf(ann, "%d %s", &i, name[offInd]);
    }

    offset[offInd] = bwt->seq_len / 2;
    char buffer[max_sam_header_size];
    buffer[0] = '\0';
    sam_headers(buffer,offset, offInd, max_sam_header_size);
    fputs(buffer, stdout);

    pthread_t *threads;
    pthread_attr_t attr;
    global_vars *data;
    threads = (pthread_t*)calloc(args->threads, sizeof(pthread_t));
    pthread_attr_init(&attr);
    data = (global_vars*)calloc(args->threads, sizeof(global_vars));
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for(j=0; j<args->threads; j++) {
        data[j].tid = j;
        data[j].bwt = bwt;
        data[j].opt = opt;
        data[j].args = args;
        data[j].ks = ks;
        data[j].ks2 = ks2;
        data[j].reference = reference;
        data[j].reference_size = reference_size;
        data[j].reference_reminder = reference_reminder;
        data[j].offset = offset;
        data[j].offInd = offInd;
    }

    running_threads_num = args->threads;
    if (args->order) {
        outBufLen = args->out_buffer_factor * args->threads * seq_num_per_read;
        outBuf = malloc(sizeof(outBuf_t) * outBufLen);
        bzero(outBuf, sizeof(outBuf_t) * outBufLen);
        outBufIndex = 1;
    }
    for (j = 0; j < args->threads; j++)
        pthread_create(&threads[j], NULL, worker2, data + j);
    if (! args->order) // Each threads writes own results in the SAM file
        for (j = 0; j < args->threads; ++j)
            pthread_join(threads[j], 0);
    else { // This thread is now writing results in the SAM file
        unsigned long long o = outBufIndex;
        while (running_threads_num > 0 || outBuf[o].read_num > 0) {
            if (!outBuf[o].read_num) usleep(sleep_time);
            else {
                if (outBuf[o].read_num != outBufIndex)
                    OutputBufferError();
                fputs(outBuf[o].str, stdout);
                free(outBuf[o].str);
                outBuf[o].read_num = 0;
                outBufIndex++;
                o = outBufIndex % outBufLen;
            }
        }
    }
    fprintf(stderr, "Alignment finished.\n");
    // destroy
    if (args->order) free(outBuf);
    free(data);
    free(threads);
    free(reference);
    free(offset);
    bwt_destroy(bwt);
    bwa_seq_close(ks);
    if (args->paired)
        bwa_seq_close(ks2);
    offInd = 0;
}

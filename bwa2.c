#include <inttypes.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include<ctype.h>
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
#include <math.h>
//#include "aligner.h"
#include "sam.h"
#define MAX(a, b) (a > b) ? a : b
#define MIN(a, b) (a < b) ? a : b
typedef struct {
    int tid;
    const gap_opt_t *opt;
    bwt_t * bwt;
    aryana_args *args;
    bwa_seqio_t * ks, *ks2;
    uint64_t * reference, reference_size, reference_reminder;
    bwtint_t * offset, offInd;
} global_vars;

struct report{
    global_vars *g;
    char * buffer;
    char ** cigar;
    char ** cigar2;
    bwa_seq_t *seq;
    bwa_seq_t *seq2;
    hash_element *table;
    hash_element * table2;
    int * can;
    int * can2;
    int canNum;
    int canNum2;
    penalty_t * penalty;
    penalty_t * penalty2;
};

void free_report(struct report* rpt){
    free(rpt->can);
    free(rpt->can2);
    free(rpt->penalty);
    free(rpt->penalty2);
    free(rpt);//it was reserved for the report struct
}

void reverse_seq(bwa_seq_t* seq){
    int jjj = seq->len - 1;
    int iii = 0;
    for (; iii < jjj && jjj >= 0 && iii < seq->len; (iii++, jjj--)) {
		ubyte_t tmp = seq->seq[iii];
		seq->seq[iii] = seq->seq[jjj];
		seq->seq[jjj] = tmp;
    }
}

int* get_seeds(aryana_args * args, int read_len, int user_seed){
    int* seeds = malloc((MAX_SEED_COUNT+1)*sizeof(int));
    if(read_len<=225){
        seeds[0] = 60, seeds[1] = 15, seeds[2] = 45;
    }else if(read_len<=275){
        seeds[0] = 60, seeds[1] = 15, seeds[2] = 45;
    }else if(read_len<=325){
        seeds[0] = 60, seeds[1] = 15, seeds[2] = 45;
    }else if(read_len<=375){
        seeds[0] = 75, seeds[1] = 15, seeds[2] = 55;
    }else if(read_len<=450){
        seeds[0] = 75, seeds[1] = 15, seeds[2] = 60;
    }else if(read_len<=550){
        seeds[0] = 100, seeds[1] = 15, seeds[2] = 60;
    }else if(read_len<=650){
        seeds[0] = 100, seeds[1] = 15, seeds[2] = 60;
    }else if(read_len<=850){
        seeds[0] = 100, seeds[1] = 15, seeds[2] = 60;
    }else{
        seeds[0] = 100, seeds[1] = 15, seeds[2] = 60; 
    }
    if(user_seed != -1){
        for(int i = MAX_SEED_COUNT; i>0; i--){
            seeds[i] = seeds[i-1];
        }
        seeds[0] = user_seed;
    }
    return seeds;
}

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


// The Data structure used for preserving the initial order of the reads
typedef struct {
    char * str;
    unsigned long long read_num;
} outBuf_t;
outBuf_t * outBuf;
int running_threads_num;
unsigned long long outBufIndex, outBufLen;

double QtoP(double quality)
{
	return pow(10,quality/(-10));
}

double PtoQ(double prob)
{
	const double eps=0.0001;
	if(prob<eps)
		return 40;
	return MIN(-10*log10(prob),40);
}

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
double prob_function(global_vars *g, bwtint_t len, ubyte_t *seq, ubyte_t * qual, char * cigar, uint64_t refrence_index)
{	
    	int read_index=0;

	const double p_snp=0.0011;	
	const double p_indel=0.0001;
	const double p_delete_open=0.003;
	const double p_delete_extend=0.05;

	const double p_insert_open=0.003;
	const double p_insert_extend=0.05;

	double p_d,p_m,p_i;

	double p_mismatch=0; // compute from read quality

	double p_ans=1;

	char prev='M';	

	int i;
    	for(i=0; cigar[i]; i++)
    	{    
   		char tmp[10];
        	int j=0;
       		while(isdigit(cigar[i]))
           		tmp[j++]=cigar[i++];
        	tmp[j]=0;
        	bwtint_t num=atoi(tmp);
		
		p_d=(prev=='D' ? p_delete_extend : p_delete_open);
		p_i=(prev=='I' ? p_insert_extend : p_insert_open);
		p_m=1-p_i-p_d;

		if(cigar[i]=='D')
		{

			prev=cigar[i];
			p_ans *= (p_d) * pow(p_delete_extend,num-1);
			refrence_index+=num;
			continue;
		}
		if(cigar[i]=='I')
		{
			p_ans *= (p_insert_open) * pow(p_insert_extend,num-1);
			prev=cigar[i];
			while(num >0)
			{
				p_ans*=QtoP(qual[read_index]-33)*p_i+(1-QtoP(qual[read_index]-33))*p_indel;
				read_index++;
				num--;
			}
			continue;
		}
		for(;num>0;read_index++, refrence_index++ ,num--)
		{

			p_mismatch=QtoP(qual[read_index]-33);
			if(seq[read_index]==getNuc(refrence_index,g->reference, g->bwt->seq_len))
			{
				p_ans *=((1-p_snp)*(1-p_mismatch)+p_snp*p_mismatch/3) * p_m;
			}
			else 
			{
				p_ans *=((1-p_snp)*(p_mismatch/3)+(p_snp/3)*(1-p_mismatch)+(2*p_snp/3)*(p_mismatch/3)) * p_m;
			}
			
			
		}
    	}
	return p_ans;
}
void compute_mapq(global_vars * g, int * can, int * can2, int num, penalty_t * p, penalty_t * p2, bwtint_t len, bwtint_t len2, ubyte_t *seq, ubyte_t *seq2, ubyte_t *qual, ubyte_t *qual2, hash_element * t, hash_element * t2, char * c[], char * c2[], int paired) {
    
    int i;
    double sum=0,sum2=0;
    for (i = 0; i < num; i++) {
        p[i].mapq = prob_function(g,len,seq, qual,c[i],t[can[i]].index);

        
	if (paired) p2[i].mapq = prob_function(g,len2, seq, qual2 ,c2[i], t2[can2[i]].index);	
	
	sum+=p[i].mapq;
	if(paired)sum2+=p2[i].mapq;
    }
    double e1=0.1;
    for (i=0; i < num; i++)
    {	
		double e2=1-p[i].mapq/sum;
		p[i].mapq=PtoQ(e1+(1-e1)*e2);//TODO prior probability
		if(paired)
		{	
			double e2=1-p2[i].mapq/sum2;
		        p2[i].mapq=PtoQ(e1+(1-e1)*e2);//TODO prior probability
		}
    }

}


void find_best_candidates(global_vars * g, int candidates_num, int candidates_num2, int candidates_size, int * candidates, int * candidates2, penalty_t * penalty, penalty_t * penalty2, int * best_candidates, int * best_candidates2, int * best_num, hash_element *table, hash_element * table2, bwtint_t len, bwtint_t len2, ubyte_t* seq, ubyte_t *seq2, ubyte_t *qual, ubyte_t *qual2, char * cigar[], char * cigar2[], int paired) { 
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
        compute_mapq(g, best_candidates, best_candidates2, *best_num, penalty, penalty2, len, len2, seq, seq2, qual, qual2, table, table2, cigar, cigar2, paired);
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
struct report* align_read(global_vars * g, char * buffer, char *cigar[], char * cigar2[], bwa_seq_t *seq, bwa_seq_t *seq2, hash_element *table, hash_element * table2, uint64_t level, int **d, char **arr, char* tmp_cigar) {
     
    buffer[0]='\0';
    int candidates_size=g->args->potents;
    int candidates[candidates_size], candidates2[candidates_size], candidates_num = candidates_size, candidates_num2 = candidates_size;
    int *best_candidates = malloc(sizeof(int)*candidates_size);
    int *best_candidates2 = malloc(sizeof(int)*candidates_size);
    int best_num;
    int best_num2;
    penalty_t* penalty = malloc(sizeof(penalty_t)*candidates_size); 
    penalty_t* penalty2 = malloc(sizeof(penalty_t)*candidates_size);

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
    find_best_candidates(g, candidates_num, candidates_num2, candidates_size, candidates, candidates2, penalty, penalty2, best_candidates, best_candidates2, &best_num, table, table2, seq->len,  (g->args->paired) ? seq2->len : 0, seq->seq, (g->args->paired) ? seq2->seq : 0, seq->qual, (g->args->paired) ? seq2->qual : 0, cigar, cigar2, g->args->paired);
    if (g->args->paired && g->args->discordant && best_num == 0) { // No valid paired alignments but we should report "discordant" alignments: the best alignment for each read of the pair separately
        find_best_candidates(g, candidates_num, 0, candidates_size, candidates, 0, penalty, 0, best_candidates, 0, &best_num, table, 0, seq->len, 0, seq->seq, seq2->seq, seq->qual, seq2->qual,  cigar, 0, 0);
        find_best_candidates(g, candidates_num2, 0, candidates_size, candidates2, 0, penalty2, 0, best_candidates2, 0, &best_num2, table2, 0, seq2->len, 0, seq->seq, seq2->seq , seq->qual, seq2->qual, cigar2, 0, 0);
        return report_discordant_alignment_results(g, buffer, cigar, cigar2, seq, seq2, table, table2, best_candidates, best_candidates2, best_num, best_num2, penalty, penalty2);
    }
    struct report* rpt = malloc(sizeof(struct report));
    rpt->g = g;
    rpt->buffer = buffer;
    rpt->cigar = cigar;
    rpt->cigar2 = cigar2;
    rpt->seq = seq;
    rpt->seq2 = seq2;
    rpt->table = table;
    rpt->table2 = table2;
    rpt->can = best_candidates;
    rpt->can2 = best_candidates2;
    rpt->canNum = best_num;
    rpt->canNum2 = best_num2;
    rpt->penalty = penalty;
    rpt->penalty2 = penalty2;
    return rpt;
}



// This is the central function of each thread. It reads a group of reads, aligns each one, and puts the output to be printed in SAM file

void multiAligner(global_vars * g) {

    // Allocating memory and initializing variables
    bwa_seq_t *seqs, *seqs2 = 0;
    int n_seqs = 0, n_seqs2;
    int total_seqs = 0;
    bwtint_t j = 0;
    hash_element *table, *table2 = 0;

    table=(hash_element *)malloc(HASH_TABLE_SIZE*(sizeof (hash_element)));
    reset_hash(table);

    int **d= malloc((MAX_READ_LEN+1)*(sizeof (int *)));
    char **arr= malloc((MAX_READ_LEN+1)*(sizeof (char *)));
    char *tmp_cigar= malloc(MAX_READ_LEN*(sizeof (char)));
    for (j=0; j<MAX_READ_LEN; j++) {
        d[j]=malloc(300*(sizeof (int)));
        arr[j]=malloc(300*(sizeof (char)));
    }
    char buffer[seq_num_per_file_read * MAX_SAM_LINE * (sizeof(char))];
    buffer[0] = '\0';
    long long lasttmpsize = 0;
    char **cigar = (char **) malloc(g->args->potents * sizeof(char *)), **cigar2 = 0;
    for (j=0; j < g->args->potents; j++)
        cigar[j]=(char *) malloc(MAX_CIGAR_SIZE*(sizeof (char)));

    if (g->args->paired) {
        table2=(hash_element *)malloc(HASH_TABLE_SIZE*(sizeof (hash_element)));
        reset_hash(table2);
        cigar2 =  (char **) malloc(g->args->potents * sizeof(char *));
        for (j=0; j < g->args->potents; j++)
            cigar2[j]=(char *) malloc(MAX_CIGAR_SIZE*(sizeof (char)));
    }
    // The thread body
    while(true) {

        pthread_mutex_lock(&input);
        int read = 1;
        if((seqs = bwa_read_seq(g->ks, seq_num_per_file_read, &n_seqs, g->opt->mode, g->opt->trim_qual)) == 0 ||
                (g->args->paired && (seqs2 = bwa_read_seq(g->ks2, seq_num_per_file_read, &n_seqs2, g->opt->mode, g->opt->trim_qual)) == 0)) read = 0;
        pthread_mutex_unlock(&input);
        if (! read) break;
        if (g->args->paired && n_seqs != n_seqs2) {
            fprintf(stderr, "Error: The number of reads from fastq1 is not equal to fastq2\n");
            return;
        }

        lasttmpsize = 0;
        int level_counter = 1;
        //------initialising variables used for choosing the best result
        int seed_count = g->args->seed_check;
        int penalties[seed_count];
        char buffs[seed_count*MAX_SAM_LINE*sizeof(char)];
        int buffs_len[seed_count];
        int user_seed = g->args->seed_length;
        for(j=0; j<n_seqs; j++) {
           
            /*reversing the sequence*/
            bwa_seq_t *seq = &seqs[j];
	        reverse_seq(seq);

            /*aligning the sequence with multiple seeds*/
            //------set the seed sets
            int *seeds = get_seeds(g->args, seq->len, user_seed);
            //------for storing the sam outputs      
            int tmpsize = 0;

            //------looping on different seeds
            int seed_counter = 0;
            for(int seed_c = 0; seed_c<seed_count; seed_c++){
                g->args->seed_length = seeds[seed_c];
                struct report* rpt = align_read(g, buffs+tmpsize, cigar, cigar2, &seqs[j], (g->args->paired)? &seqs2[j] : 0, table, table2,total_seqs + level_counter, d,arr, tmp_cigar);
                level_counter++;
                penalties[seed_c] = (rpt->canNum?rpt->penalty[0].penalty:maxPenalty);
                seed_counter++;
                int tmp=report_alignment_results(rpt->g, rpt->buffer, rpt->cigar,rpt->cigar2, rpt->seq,rpt->seq2, rpt->table, rpt->table2, rpt->can, rpt->can2, rpt->canNum, rpt->canNum2, rpt->penalty, rpt->penalty2);
                buffs_len[seed_c] = tmp;
                tmpsize+=tmp;
                /*free the memory reserved for the aligner*/
                free_report(rpt);
                /*if it is a good one, escape so*/
                if(penalties[seed_c] < g->args->gap_ext_penalty*0.05*seq->len)
                    break;
            }
            //------free the seeds array
            free(seeds);
            
            //------choose the best alignment
            int minPenalty = maxPenalty, minIndex = 0;
            for(int seed_c = 0; seed_c<seed_counter; seed_c++){
                if(penalties[seed_c] < minPenalty){
                    minPenalty = penalties[seed_c];
                    minIndex = seed_c;
                }
            }
            //------copy the best alignment sam file result to buffer
            int where = 0;
            for(int seed_c = 0; seed_c<minIndex; seed_c++){
                where+=buffs_len[seed_c];
            }
            memcpy(buffer+lasttmpsize, buffs+where, buffs_len[minIndex]);
            lasttmpsize+=buffs_len[minIndex];
            
           

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
                if (j % seq_num_per_file_read == 0 || j == n_seqs - 1) {
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
    for (j=0; j<MAX_READ_LEN; j++) {
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
    fprintf(stderr, "Thread %d finished.\n", g->tid);
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

void byebye(char * message) {
	fprintf(stderr, "%s\n", message);
	exit(1);
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
    if (! fn_fa) {
        fprintf(stderr, "No input file is provided\n");
        exit(1);
    }
    ks = bwa_open_reads(opt->mode, fn_fa);
    if (args->paired) {
        if (! fn_fa2) {
            fprintf(stderr, "Paired input file is not provided\n");
            exit(1);
        }
        ks2 = bwa_open_reads(opt->mode, fn_fa2);
    }
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
    if (!fgets(line, sizeof line, ann) || !fscanf(ann, "%d %s", &i, name[0])) {
        fclose(ann);
		bwa_seq_close(ks);
		if (args->paired)
        	bwa_seq_close(ks2);
		byebye("Error reading one of the index files");
    }
	while(fgets(line, sizeof line, ann) != NULL) {
        if (! fscanf(ann, "%llu", (unsigned long long *) &offset[offInd++]) ||
			!fgets(line, sizeof line, ann) || 
			!fscanf(ann, "%d %s", &i, name[offInd])) {
            fclose(ann);
            bwa_seq_close(ks);
        	if (args->paired)
            	bwa_seq_close(ks2);
			byebye("Error reading one of the index files");
		}
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
        outBufLen = args->out_buffer_factor * args->threads * seq_num_per_file_read;
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
    fclose(ann);
    bwt_destroy(bwt);
    bwa_seq_close(ks);
    if (args->paired)
        bwa_seq_close(ks2);
    offInd = 0;
}


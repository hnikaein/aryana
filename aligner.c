#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>

//#include "bwtgap.h"
//#include "bwtaln.h"
#include "aligner.h"
//#include "utils.h"
#include "smith.h"
#include <math.h>
#define false 0
#define true 1
#define bool int
#include "bwt2.h"
extern int debug;
extern void PrintSeq(const unsigned char *, int, int);
extern void PrintRefSeq(unsigned long long, unsigned long long, unsigned long long, int);
extern void GetRefSeq(unsigned long long, unsigned long long, unsigned long long, char*);
const char atom[4]= {'A','C','G','T'};
extern char getNuc(uint64_t place, uint64_t seq_len);

int create_cigar(hash_element *best, char *cigar, int len, const ubyte_t *seq, uint64_t seq_len,int **d, char **arr, char * tmp_cigar )
{
    int *valid=(int *)malloc(best->parts*(sizeof (int)));
    bwtint_t lastvalid=best->parts-1;
    long long i=0,j=0;
    for (i=best->parts-1; i>=0; i--)
    {
        valid[i]=1;
        if (abs((best->match_index[i]-best->match_start[i])-best->index) > 5+best->match_start[i]/20)
            valid[i]=0;
        if (valid[i]==1 && i<best->parts-1 && (best->match_start[lastvalid]+best->matched[lastvalid]>best->match_start[i]))
        {
            valid[i]=0;
            if (abs((best->match_index[i]-best->match_start[i])-best->index) < abs((best->match_index[lastvalid]-best->match_start[lastvalid])-best->index))
            {
                valid[lastvalid]=0;
                valid[i]=1;
            }
        }
        if (valid[i]==1)
            lastvalid=i;
        if (i==0)
            break;
    }

    //for (i=0; i<best->parts; i++)
    //	fprintf(stderr, "part %lld: %lld %lld %d %d\n",i,best->match_start[i],best->match_index[i],best->matched[i],valid[i]);

    int slack=20;
    bwtint_t head_match=0,head_index=best->index >= slack ? best->index-slack : 0;
    bwtint_t slack_index=head_index;
    int print_head=0;
    int total_errors=0;
    if (best->parts<=0 || best->parts > 50)
        fprintf(stderr , "too much parts!\n");
    if (debug > 0) {
        fprintf(stderr, "Generating Cigar, best->parts:%d, Seqs:\n", best->parts);
        PrintSeq(seq, len, 1);
        PrintRefSeq(head_index, head_index+len+2*slack, seq_len, 1);
    }
    if (0)
    {
        int errors_=0;
        print_head=smith_waterman(0, len, head_index, head_index+len+2*slack, cigar, print_head, seq, len, &errors_, seq_len, d, arr, tmp_cigar);
        total_errors+=errors_;
        if (debug > 1) {
            fprintf(stderr, "SmithWaterman1 [%d,%d],[%llu,%llu] errors %d, cigar %s, tmp_cigar %s\n", 0, len, (unsigned long long) head_index, (unsigned long long) head_index + len + 2*slack, errors_, cigar, tmp_cigar);
            PrintSeq(seq, len, 1);
            PrintRefSeq(head_index, head_index+len+2*slack, seq_len, 1);
        }

    }
    else
        for (i=best->parts-1; i>=0; i--)
        {
            if (valid[i] && !(best->match_start[i] < head_match || best->match_index[i] < head_index))// && abs((best->match_index[i]-head_index)-(best->match_start[i]-head_match)<=3+(best->match_index[i]-head_index)/20))
            {
                int errors=0;
                print_head=smith_waterman(head_match, best->match_start[i], head_index, best->match_index[i], cigar, print_head, seq, len, &errors, seq_len, d, arr, tmp_cigar);
                //fprintf(stderr,"start: %llu, end: %llu, errors :: %d\n",head_match,best->match_start[i],errors);
                total_errors+=errors;
                if (debug > 1) {
                    fprintf(stderr, "SmithWaterman2 [%llu,%llu],[%llu,%llu] errors %d, cigar %s, tmp_cigar %s\n", (unsigned long long) head_match, (unsigned long long) best->match_start[i], (unsigned long long) head_index, (unsigned long long) best->match_index[i], errors, cigar, tmp_cigar);
                    PrintSeq(seq+head_match,  best->match_start[i] - head_match + 1, 1);
                    PrintRefSeq(head_index, best->match_index[i], seq_len, 1);
                }
                print_head+=snprintf(cigar+print_head,10,"%"PRIu64"%c",best->matched[i],'m');
                head_match=best->match_start[i]+best->matched[i];
                head_index=best->match_index[i]+best->matched[i];
            }
            if (i==0)
            {
                int errors=0;
                print_head=smith_waterman(head_match,len,head_index, head_index+len-head_match+slack, cigar, print_head,seq, len, &errors, seq_len, d, arr, tmp_cigar);
                //fprintf(stderr,"start: %llu, end: %llu, errors :: %d\n",head_match,len,errors);
                total_errors+=errors;
                if (debug > 1) {
                    fprintf(stderr, "SmithWaterman3 [%llu,%d],[%llu,%llu] errors %d, cigar %s, tmp_cigar %s\n", (unsigned long long) head_match, len, (unsigned long long) head_index, (unsigned long long) head_index+len-head_match+slack, errors, cigar, tmp_cigar);
                    PrintSeq(seq+head_match, len - head_match + 1, 1);
                    PrintRefSeq(head_index, head_index+len-head_match+slack, seq_len, 1);
                }
                break;
            }
        }

    //refining cigar
    int firstblood=1;
    int last_size=0;
    char last_char='m';
    print_head=0;
    for (i=0; cigar[i]; i++)
    {
        char tmp[10];
        j=0;
        while(isdigit(cigar[i]))
            tmp[j++]=cigar[i++];
        tmp[j]=0;
        bwtint_t num=atoi(tmp);
        if (firstblood==1)
        {
            firstblood=0;

            last_size=num;
            last_char=cigar[i];

            if (cigar[i]=='d')
            {
                last_size=0;
                slack_index+=num;
            }

            continue;
        }
        if (cigar[i]==last_char)
        {
            if (last_size>0)
                last_size+=num;
        }
        else
        {
            if (last_size>0)
                print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);
            last_size=num;
            last_char=cigar[i];
        }
    }
    if (last_char!='d')
        print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);

    best->index=slack_index;
    free(valid);
    //fprintf(stderr,"total errors:: %d\n",total_errors);
    return total_errors;//best->index;
}

#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)

void floyd(long long down, long long up , int exactmatch_num,long long *selected){
    long N = up - down+1;
    long long in, im,j;
    im = 0;
    if(N < exactmatch_num){
        for (j=down; j<=up; j++) {
            selected[j-down]=j;
        }
        return;
    }
    char is_used[BITNSLOTS(N)];
    memset(is_used, 0, BITNSLOTS(N));

    in = N - exactmatch_num;
    
    for (; in < N && im < exactmatch_num; ++in) {
        srand (time(NULL));
        long long r = rand() % (in + 1); /* generate a random number 'r' */
        if(BITTEST(is_used, r))
        /* we already have 'r' */
            r = in; /* use 'in' instead of the generated number */

        selected[im++] = r + down; /* +1 since your range begins from 1 */
        BITSET(is_used, r);
        
    }
}
void knuth(long long down, long long up , int exactmatch_num,long long *selected){
    long long j,in, im=0 ,rn,rm;
    long N = up - down+1;
    if(N < exactmatch_num){
        for (j=down; j<=up; j++) {
            selected[j-down]=j;
        }
        return;
    }
    for (in = 0; in < N && im < exactmatch_num; ++in) {
        rn = N - in;
        rm = exactmatch_num - im;
        srand (time(NULL));
        if (rand() % rn < rm)
            selected[im++] = down + in;
    }
}
void match_select(long long down, long long up , int exactmatch_num,long long *selected){
    if (exactmatch_num < (up-down+1)/2) {
        floyd(down,up,exactmatch_num, selected);
    }
    else
        knuth(down,up,exactmatch_num, selected);
}

// The main Aryana aligner routine.

void aligner(bwt_t *const bwt, int len, ubyte_t *seq, bwtint_t level, hash_element * table, int *best, int best_size, int *best_found, aryana_args *args)
{
    if (debug > 4) {
        fprintf(stderr, "Generating BWT Table, seq_len= %lld...", (long long) bwt->seq_len);
        long long l, i;
        FILE * f = fopen("BWT.txt", "w");
        for (l = 0; l < bwt->seq_len; l++) {
            i = bwt_sa(bwt, l);
            fprintf(f, "%lld\n", i);
        }
        fclose(f);
        fprintf(stderr, "Done.\n");
        exit(0);
    }
    //showerr(len,seq);
    //initialize
//	fprintf(stderr,"paired :: %d\n", options.paired);
    bwtint_t down, up;
    bwtint_t limit;
    long long i = 0, j = 0;

    /*best->index = best->value = best->place = best->level = 0;
    best->parts = 0;
    best->index = bwt->seq_len;
    second->index = second->value = second->place = second->level = 0;
    second->index = bwt->seq_len;
    second->parts=0;
    */
    //for (i=0; i<100; i++)
    //	fprintf(stderr,"%c",atom[getNuc(10000+i,bwt->seq_len)]);
    //fprintf(stderr,"\n");
    //set k
    bwtint_t k;
    if (args->seed_length==-1)
    {
        k = 26;
        if (len < 600)
            k=24;
        if (len < 300)
            k=22;
        if (len < 150)
            k=20;
        if (len <80)
            k=18;
        if (len <40)
            k=15;
    }
    else
        k=args->seed_length;

    //reversing
    j = len-1;
    i=0;
    for (; i<j && j >= 0 && i < len; (i++,j--)) {
        ubyte_t tmp = seq[i];
        seq[i] = seq[j];
        seq[j] = tmp;
        //if (seq[i]>3 || seq[j]>3)
        //	return;
    }
    //inexact match
    bwtint_t groupid_last=1;
    for(i=len - 1; i>=k; i--) {
        bwt_match_limit_rev(bwt, k, seq+i - k + 1, &down, &up,&limit);
        if(limit < k) {
            i = i - k + limit;
            if (i < k) break;
            continue;
        }
        bwt_match_limit(bwt, i+1, seq, &down, &up,&limit);
        if (debug > 2) {
            fprintf(stderr, "aligner(), %llu regions have exact match with %llu score, seq: ", (unsigned long long) (up - down + 1), (unsigned long long) limit);
            PrintSeq(seq + i + 1 - limit, limit, 1);
        }

        long long selected[args->exactmatch_num];
        match_select(down, up , args->exactmatch_num,selected);
        //		for (j=down; j<=up && j <= (down + 50); j++){
        int ii=0;
        for ( j= selected[ii]; ii < args->exactmatch_num && ii < up ; ii++){
            bwtint_t index=bwt_sa(bwt,j);
//            		if (index > len)
//                 		index = index - len;
//             		else
//                 		index = 0;
            if (debug > 2) {
                fprintf(stderr, "match, index: %llu, seq: ", (unsigned long long) index);
                PrintRefSeq(index, index + limit - 1, bwt->seq_len, 1);
                if (debug > 3) {
                    fprintf(stderr, "Additional sequences:\n");
                    long long l;
                    for (l = j - 3; l < j + 3; l++) {
                        if (l < 0 || l > bwt->seq_len) continue;
                        bwtint_t tmpind = bwt_sa(bwt, l);
                        PrintRefSeq(tmpind, tmpind + limit - 1, bwt->seq_len, 1);
                    }
                }
            }
            bwtint_t score=limit;

            bwtint_t rindex=index- (i - limit+1);
            add(bwt, rindex/len, score, level, index - (i - limit+1), best, best_size, best_found, table, i - limit+1, limit, index,len, groupid_last); // if level changed, check the find_value in hash.c
        }
        groupid_last++;
        if(i > k) {
            if((limit - k + 1) > 0)
                i = i - limit + (k - 1);
            else
                fprintf(stderr, "manfi\n");
            if (i < k) break;
        }
    }
    if (args->best_factor!=-1)
    {
        for (i=best_size-2; i>=(*best_found); i--)
        {
            if (best[i] == -1)
                break;
            if (best[i]<best[best_size-1]/args->best_factor)
            {
                (*best_found)=i+1;
                break;
            }
            if (i==0)
                break;
        }
    }

}


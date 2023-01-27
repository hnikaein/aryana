#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <assert.h>
#include "hash.h"
#include "aryana_args.h"
#include "bwa2.h"
#include "smith.h"
//#include "bwtgap.h"
//#include "bwtaln.h"
#include "aligner.h"
//#include "utils.h"
#include <math.h>
#define false 0
#define true 1
#define bool int
#include "bwt2.h"
extern int debug;
extern void PrintSeq(const unsigned char *, int, int);
extern void PrintRefSeq(uint64_t * reference, unsigned long long, unsigned long long, unsigned long long, int);
extern void GetRefSeq(unsigned long long, unsigned long long, unsigned long long, char*);
const char atom[4]= {'A','C','G','T'};
extern char getNuc(uint64_t place, uint64_t * reference, uint64_t seq_len);

extern long long total_candidates, best_factor_candidates;

void create_cigar(aryana_args * args, hash_element *best, char *cigar, int len, const ubyte_t *seq, uint64_t seq_len,int **d, char **arr, char * tmp_cigar, penalty_t * penalty, uint64_t * reference, ignore_mismatch_t ignore)
{
    penalty->mismatch_num = 0;
    penalty->gap_open_num = 0;
    penalty->gap_ext_num = 0;
    int *valid=(int *)malloc(best->parts*(sizeof (int)));
    bwtint_t lastvalid=best->parts-1;
    long long i=0,j=0;
    for (i=best->parts-1; i>=0; i--)
    {
        valid[i]=1;
        if (llabs((signed)(best->match_index[i]-best->match_start[i])-(signed)best->index) > 5+best->match_start[i]/20)
            valid[i]=0;
        if (valid[i]==1 && i<best->parts-1 && (best->match_start[lastvalid]+best->matched[lastvalid]>best->match_start[i]))
        {
            valid[i]=0;
            if (llabs((signed)(best->match_index[i]-best->match_start[i])-(signed)best->index) < llabs((signed)(best->match_index[lastvalid]-best->match_start[lastvalid])-(signed)best->index))
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

    int slack=10;
    bwtint_t head_match=0,head_index=best->index >= slack ? best->index-slack : 0;
    bwtint_t slack_index=head_index;
    int print_head=0;
    if (best->parts<=0 || best->parts > 50)
        fprintf(stderr , "too much parts!\n");
    if (args->debug > 0) {

        fprintf(stderr, "Generating Cigar, best->parts:%d, Seqs:\n", best->parts);
        PrintSeq(seq, len, 1);
        PrintRefSeq(reference, head_index, head_index+len+2*slack, seq_len, 1);
    }
    for (i=best->parts-1; i>=0; i--)
    {
        if (valid[i] && !(best->match_start[i] < head_match || best->match_index[i] < head_index))// && llabs((signed)(best->match_index[i]-head_index)-(signed)(best->match_start[i]-head_match)<=3+(best->match_index[i]-head_index)/20))
        {
            print_head=smith_waterman(args, head_match, best->match_start[i], head_index, best->match_index[i], cigar, print_head, seq, len, &penalty->mismatch_num, seq_len, d, arr, tmp_cigar, reference, ignore);
            if (args->debug > 1) {
                fprintf(stderr, "SmithWaterman1 [%llu,%llu],[%llu,%llu] cigar %s, tmp_cigar %s\n", (unsigned long long) head_match, (unsigned long long) best->match_start[i], (unsigned long long) head_index, (unsigned long long) best->match_index[i], cigar, tmp_cigar);
                PrintSeq(seq+head_match,  best->match_start[i] - head_match, 1);
                PrintRefSeq(reference, head_index, best->match_index[i]-1, seq_len, 1);
            }
            print_head+=snprintf(cigar+print_head,10,"%"PRIu64"%c",best->matched[i],'M');
            head_match=best->match_start[i]+best->matched[i];
            head_index=best->match_index[i]+best->matched[i];
        }
        if (i==0)
        {
            bwtint_t end = head_index+len-head_match+slack;
            if (end >= seq_len) end = seq_len -1;
            if (head_index > end || (signed) (len - head_match) > (signed) (end - head_index) + 10)
                print_head+=snprintf(cigar+print_head,10,"%"PRIu64"%c",(len-head_match),'I');
            else {
                print_head=smith_waterman(args, head_match,len,head_index, end, cigar, print_head,seq, len, &penalty->mismatch_num, seq_len, d, arr, tmp_cigar, reference, ignore);
                //fprintf(stderr,"start: %llu, end: %llu, errors :: %d\n",head_match,len,errors);
                if (args->debug > 1) {
                    fprintf(stderr, "SmithWaterman2 [%llu,%llu],[%llu,%llu] cigar %s, tmp_cigar %s\n", (unsigned long long) head_match, (unsigned long long) len, (unsigned long long) head_index, (unsigned long long) end, cigar, tmp_cigar);
                    PrintSeq(seq+head_match, len - head_match, 1);
                    PrintRefSeq(reference, head_index, end-1, seq_len, 1);
                }
            }
            break;
        }
    }

    //refining cigar
    int firstblood=1;
    int last_size=0;
    char last_char='M';
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
            last_size=num;
            last_char=cigar[i];
            if (cigar[i]=='D')
            {
                last_size=0;
                slack_index+=num;
            } else firstblood = 0;
            continue;
        }
        if (cigar[i]==last_char)
        {
            if (last_size>0)
                last_size+=num;
        }
        else
        {
            if (last_size>0) {
                print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);
                if (last_char == 'I' || last_char == 'D') {
                    penalty->gap_open_num++;
                    penalty->gap_ext_num += last_size - 1;
                }
            }
            last_size=num;
            last_char=cigar[i];
        }
    }
    if (last_char!='D') {
        if (last_size > 0)
            print_head+=snprintf(cigar+print_head,10,"%d%c",last_size,last_char);
        if (last_char == 'I') {
            penalty->gap_open_num++;
            penalty->gap_ext_num += last_size - 1;
        }
    }
    if (args->debug > 0)
        fprintf(stderr, "Cigar: %s, Mismatch: %d, Gap Open: %d, Gap Ext: %d\n", cigar, penalty->mismatch_num, penalty->gap_open_num, penalty->gap_ext_num);
    best->index=slack_index;
    free(valid);
}

#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)

void floyd(long long down, long long up , int exactmatch_num,long long *selected) {
    long N = up - down+1;
    long long in, im,j;
    im = 0;
    if(N < exactmatch_num) {
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

void knuth(long long down, long long up , int exactmatch_num,long long *selected) {
    long long j,in, im=0 ,rn,rm;
    long N = up - down+1;
    if(N < exactmatch_num) {
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
void match_select(long long down, long long up , int exactmatch_num,long long *selected) {
    if (exactmatch_num < (up-down+1)/2) {
        floyd(down,up,exactmatch_num, selected);
    }
    else
        knuth(down,up,exactmatch_num, selected);
}

// The main Aryana aligner routine.

void aligner(bwt_t *const bwt, int len, ubyte_t *seq, bwtint_t level, hash_element * table, int *best, int best_size, int *best_found, aryana_args *args, uint64_t * reference)
{
    // initialize
    bwtint_t down, up;
    bwtint_t limit;
    long long i = 0, j = 0;

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
        if (args->bisulfite)
            k = k * 2 / 3;
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
            i = i - k + limit + 1;
            continue;
        }
        bwt_match_limit(bwt, i+1, seq, &down, &up,&limit);
        if (args->debug > 2) {
            fprintf(stderr, "aligner(), %llu regions have exact match with %llu score, seq: ", (unsigned long long) (up - down + 1), (unsigned long long) limit);
            PrintSeq(seq + i + 1 - limit, limit, 1);
        }

        long long selected[args->exactmatch_num];
        match_select(down, up , args->exactmatch_num,selected);
        //		for (j=down; j<=up && j <= (down + 50); j++){
        int ii=0;
        for ( ; ii < args->exactmatch_num && ii < up-down+1 ; ii++) {
            j= selected[ii];
            bwtint_t index=bwt_sa(bwt,j);
//            		if (index > len)
//                 		index = index - len;
//             		else
//                 		index = 0;
            if (args->debug > 2) {
                fprintf(stderr, "match, index: %llu, seq: ", (unsigned long long) index);
                PrintRefSeq(reference, index, index + limit - 1, bwt->seq_len, 1);
                if (args->debug > 3) {
                    fprintf(stderr, "Additional sequences:\n");
                    long long l;
                    for (l = j - 3; l < j + 3; l++) {
                        if (l < 0 || l > bwt->seq_len) continue;
                        bwtint_t tmpind = bwt_sa(bwt, l);
                        PrintRefSeq(reference, tmpind, tmpind + limit - 1, bwt->seq_len, 1);
                    }
                }
            }
            bwtint_t score=limit;
            if (index < (i-limit+1)) continue;
            bwtint_t rindex=index- (i - limit+1);
            assert(rindex < bwt->seq_len);
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
	total_candidates += best_size - *best_found;
    if (args->best_factor > 0)
    {
        for (i=best_size-2; i>=(*best_found); i--)
        {
            if (best[i] == -1)
                break;
            if (table[best[i]].value<table[best[best_size-1]].value*args->best_factor)
            {
				best_factor_candidates += i + 1 - *best_found;
                (*best_found)=i+1;
                break;
            }
            if (i==0)
                break;
        }
    }

}


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aryana_args.h"
#include "bwa2.h"
#include "smith.h"
int max(int q , int p)
{
    return q<p ? p : q;
}


// Adds one letter *edit* to the cigar sequence *tmp_cigar* and updates the mismatch/gap open/gap extension numbers

void extend_cigar(char * tmp_cigar, char edit, char *last, int * tail, penalty_t * penalty) {
	switch (edit) {
		case 'm': break;
		case 'M': edit = 'm'; penalty->mismatch_num++; break;
		case 'd': case 'i': if (edit != *last) penalty->gap_open_num++; else penalty->gap_ext_num++; break;
	}
	*last = edit;
	tmp_cigar[(*tail)++] = edit;
}


// Just increments the mismatch_num, gap_open_num and gap_ext_num. They might have previous values and those values are not cleared.
int smith_waterman(uint64_t match_start, uint64_t match_end, uint64_t index_start, uint64_t index_end, char *cigar, int head, const ubyte_t *read, int len, penalty_t * penalty, uint64_t seq_len, int **d, char **arr, char *tmp_cigar, uint64_t * reference)
{
//	fprintf(stderr,"read: %llu - %llu, ref: %llu - %llu\n",match_start,match_end,index_start,index_end);
    int off=max(2*(abs((index_end-index_start)-(match_end-match_start))),10);
    if (off > 98)
    {
        fprintf(stderr,"off too long\n");
        head+=snprintf(cigar+head,10,"%d",abs((index_end-index_start)-(match_end-match_start)));
        if (match_end-match_start > index_end-index_start)
        {
            head+=snprintf(cigar+head,10,"i");
            head+=snprintf(cigar+head,10,"%"PRIu64"m",index_end-index_start);
        }
        else
        {
            head+=snprintf(cigar+head,10,"d");
            head+=snprintf(cigar+head,10,"%"PRIu64"m",match_end-match_start);
        }
        return head;
    }
    off*=2;
    if (off>98)
        off=100;
//	off=200;


    long long i=0,j=0;
    for (i=off/2; i<off; i++)
    {
        if (match_start==0 || match_start==len)
            d[0][i]=0;
        else
            d[0][i]=i-off/2;
    }

    if (index_end-index_start==0 && match_end-match_start==0)
        return head;

    //char atomic[4] = { 'A' , 'C' , 'G' , 'T'};
    for (i=1; i<=match_end-match_start; i++)
    {
        for (j=0; j<off; j++)
        {
            /*			if (i<= 0 || i>6000)
            				fprintf(stderr,"i :: %llu\n",i);
            			if (j <0 || j>100)
            				fprintf(stderr,"j :: %llu\n",i);
            */			int real_off=j-off/2;
            if (real_off<0 && i < -real_off )
                continue;
            uint64_t ref_i=i+real_off;
            if (ref_i==0)
            {
                d[i][j]=i;
                continue;
            }

            d[i][j]=d[i-1][j];
            arr[i][j]='m'; // match

            if (index_start+ref_i-1 >= seq_len) { 
				d[i][j]++;
				arr[i][j] = 'd';
			}

			if (getNuc(index_start+ref_i-1,reference, seq_len)!=read[match_start+i-1]){
                d[i][j]++;
				arr[i][j] = 'M'; // mismatch
			}

            if (match_end==len && i==match_end-match_start)
            {
                if (j>0 && (d[i][j] > d[i][j-1]) )
                {
                    d[i][j]=d[i][j-1];
                    arr[i][j]='d';
                }
            }
            else
            {
                if (j>0 && (d[i][j] > d[i][j-1]+1) )
                {
                    d[i][j]=d[i][j-1]+1;
                    arr[i][j]='d';
                }
            }
            if (j<off-1 && (d[i][j] > d[i-1][j+1]+1))
            {
                d[i][j]=d[i-1][j+1]+1;
                arr[i][j]='i';
            }
        }
    }
    int cur_off=(index_end-index_start)-(match_end-match_start)+off/2;
    int cur_i=match_end-match_start;
    int tail=0;
	char last = 'm';
    while(1)
    {
        int ref_i=cur_i+cur_off-off/2;
        if (cur_i==0)
        {
            for (j=0; j<ref_i; j++)
				extend_cigar(tmp_cigar, 'd', &last, &tail, penalty);
            break;
        }
        if (ref_i==0)
        {

            for (j=0; j<cur_i; j++)
				 extend_cigar(tmp_cigar, 'i', &last, &tail, penalty);
            break;
        }

        extend_cigar(tmp_cigar, arr[cur_i][cur_off], &last, &tail, penalty);
        if (arr[cur_i][cur_off]=='i')
        {
            cur_off++;
            cur_i--;
        }
        else if (arr[cur_i][cur_off]=='d')
        {
            cur_off--;
        }
        else if (arr[cur_i][cur_off]=='m' || arr[cur_i][cur_off]=='M')
        {
            cur_i--;
        }
        else
        {
            fprintf(stderr,"There is something wrong\n");
            exit(0);
        }
    }
    tmp_cigar[tail]=0;
    int ct=0;
    last=tmp_cigar[tail-1];
    for (i=tail-1; i>=0; i--)
    {
        if (tmp_cigar[i]==last )
            ct++;
        else
        {
            head += snprintf(cigar+head,10,"%d%c",ct,last);
            last=tmp_cigar[i];
            ct=1;
        }
        if (i==0)
            break;
    }
    head+=snprintf(cigar+head,10,"%d%c",ct,last);
//	fprintf(stderr, "cigar at smith: %s\n",cigar);
    return head;
}


#include <stdio.h>
#include <inttypes.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <stdint.h>
#include <pthread.h>

#include "sam.h"
const int sam_line=5000;
char int_to_bp[4]= {'A','C','G','T'};
static char star[2] = "*";
static char equal[2] = "=";

uint64_t reverse_cigar(char *cigar)
{
    int i=0;
    int cigar_size=0;
        cigar_size++;
    char *cigar_tmp=(char *)malloc((cigar_size+50)*(sizeof (char)));
    int num=0;
    char type=cigar[cigar_size-1];
    int offset=1;
    int lasttmpsize=0;
    int length=0;
    for (i=cigar_size-2; i>=0; i--)
    {
        if (cigar[i]>='0' && cigar[i]<='9' )
        {
            num+=offset*(cigar[i]-'0');
            offset*=10;
        }
        else
        {
            lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
            if (type=='m' || type=='d')
                length+=num;
            type=cigar[i];
            num=0;
            offset=1;
        }
    }
    lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
    if (type=='m')
        length+=num;
    for (i=0; i<cigar_size; i++)
        cigar[i]=cigar_tmp[i];
    free(cigar_tmp);
    return length;
}

int binary_search(int bot, int top, uint64_t * a, uint64_t key) {
    while (bot <= top) {
        int mid = (bot + top) / 2;
        //fprintf(stderr, "bot = %d, top = %d, a[%d] = %llu\n", bot, top, mid, a[mid]);
        if(key >= a[mid] && (mid + 1 > top || a[mid + 1] > key))
            return mid;
        if(key < a[mid]) top = mid - 1;
        else bot = mid + 1;
    }
    return -1;
}

void get_genomic_position(int * flag, bwtint_t index, genomic_location * loc, char ** cigar, int * len, uint32_t * mapq, bwtint_t seq_len, int unmapped, int reverse, bwtint_t * offset, bwtint_t offInd) {
    if (index <= 0 || index >= seq_len)
        *flag |= unmapped;
    if((*flag) & unmapped) {
        mapq = 0;
        loc->pos = 0;
        loc->rname = star;
        *cigar = star;
        *len = 0;
    } else {
        if (index >= seq_len / 2) {
            index = (seq_len / 2) - (index - (seq_len / 2))-reverse_cigar(*cigar);
            *flag |= reverse;
        }
		bwtint_t chr =  binary_search(0, offInd - 1, offset, index);
		index -= offset[chr];
		loc->pos = index;
		loc->rname = name[chr];
    }       
}

int sam_generator(char *buffer, char *qname, int flag, uint32_t mapq, bwtint_t index, char * cigar, bwtint_t index2, char * cigar2, ubyte_t *seq, ubyte_t * quality, int len, int len2, bwt_t * bwt, bwtint_t * offset, bwtint_t offInd)
{
	if (! cigar) cigar = star;
	if (! cigar2) cigar2 = star;
	genomic_location first, next;
	long long tlen = 0; // Template length: the size of fragment for paired-end reads
    char  seq_string[len+1], quality_string[len+1];
	uint32_t tmp;
	get_genomic_position(&flag, index, &first, &cigar, &len, &mapq, bwt->seq_len, sf_unmapped, sf_reverse, offset, offInd);
	if (flag & sf_paired) get_genomic_position(&flag, index2, &next, &cigar2, &len2, &tmp, bwt->seq_len, sf_mate_unmapped, sf_mate_reverse, offset, offInd);
	else {
		next.pos = 0; 
		next.rname = star;
	}
	if (flag & sf_pair_mapped) {
		if (first.pos < next.pos) tlen = next.pos + len2 - first.pos;
		else tlen = -(first.pos + len - next.pos);
	}
	if ((flag & sf_paired) && ((flag & (sf_unmapped | sf_mate_unmapped)) == 0) && (strcmp(first.rname, next.rname)!=0) && strcmp(first.rname, "*")) next.rname = equal;

	seq_string[len]=quality_string[len]=0;
    int i;
    if (! (flag & sf_reverse)) {
        for (i=0; i<len; i++)
        {
            if (seq[i]<4)
                seq_string[i]=int_to_bp[seq[i]];
            else
                seq_string[i]='N';
            quality_string[i]=quality[i];
        }
    } else {
        for (i=0; i<len; i++)
        {
            if (seq[i]<4)
                seq_string[len - i - 1]=int_to_bp[3 - seq[i]];
            else
                seq_string[len - i - 1]='N';
            quality_string[len - i - 1]=quality[i];
        }
    }
    buffer[0] = '\0';
    return snprintf(buffer,sam_line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%"PRIu64"\t%lld\t%s\t%s\n",qname,flag,first.rname, first.pos, mapq, cigar,next.rname,next.pos,tlen,seq_string,quality_string);
}

int sam_headers(char * buffer, bwtint_t *  offset, int size) {
    int head=0;
    int i;
    //fprintf(stderr, "offset[size] = %llu, offset[size - 1] = %llu, last = %llu, size = %d\n", offset[size], offset[size - 1], offset[size] - offset[size - 1], size);
    head+=snprintf(buffer,40,"@HD\tVN:1.0\tSO:unsorted\n");
    //fprintf(stderr, "buffer = %s offset[0] = %llu\n name[0] = %s\n", buffer, offset[0], name[0]);
    for(i=0; i<size; i++) { //TODO
        head+=snprintf(buffer+head,40,"@SQ\tSN:%s\tLN:%"PRIu64"\n",name[i], offset[i + 1] - offset[i]); //TODO
    }
    return head;
    //return buffer;
}

/*int main()
{
	char buf[400];
	//int head=0;
	char *names[2]={"ref1","ref2"};
	long long offset[2]={25,44};
	int head=sam_headers(buf,names,offset,2);
	sam_generator(buf,head,"qname",16,"ref",53,24,"cigar","rnext",999,19,"ACTF","#FDF#D");
	printf("%s",buf);
	return 0;
}
*/


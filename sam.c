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
#include "const.h"
#include "sam.h"

char int_to_bp[4] = {'A', 'C', 'G', 'T'};
static char star[2] = "*";
static char equal[2] = "=";

uint64_t cigar_len(char *cigar) {
    int i = 0, cigar_size = strlen(cigar);
    int num = 0;
    char type = cigar[cigar_size - 1];
    int offset = 1;
    int length = 0;
    for (i = cigar_size - 2; i >= 0; i--) {
        if (cigar[i] >= '0' && cigar[i] <= '9') {
            num += offset * (cigar[i] - '0');
            offset *= 10;
        } else {
            if (type == 'M' || type == 'D')
                length += num;
            type = cigar[i];
            num = 0;
            offset = 1;
        }
    }
    if (type == 'M')
        length += num;
    return length;
}

void reverse_cigar(char *cigar) {
    int i = 0, cigar_size = strlen(cigar);
    char cigar_tmp[cigar_size + 50];
    int num = 0;
    char type = cigar[cigar_size - 1];
    int offset = 1;
    int lasttmpsize = 0;
    for (i = cigar_size - 2; i >= 0; i--) {
        if (cigar[i] >= '0' && cigar[i] <= '9') {
            num += offset * (cigar[i] - '0');
            offset *= 10;
        } else {
            lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num, type);
            type = cigar[i];
            num = 0;
            offset = 1;
        }
    }
    lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num, type);
    memcpy(cigar, cigar_tmp, cigar_size);
}

int binary_search(int bot, int top, uint64_t *a, uint64_t key) {
    while (bot <= top) {
        int mid = (bot + top) / 2;
        //fprintf(stderr, "bot = %d, top = %d, a[%d] = %llu\n", bot, top, mid, a[mid]);
        if (key >= a[mid] && (mid + 1 > top || a[mid + 1] > key))
            return mid;
        if (key < a[mid]) top = mid - 1;
        else bot = mid + 1;
    }
    return -1;
}

int get_chr(bwtint_t index, bwtint_t seq_len, bwtint_t *offset, bwtint_t offInd) {
    if (index <= 0 || index >= seq_len)
        return -1;
    return binary_search(0, offInd - 1, offset, index);
}

void
get_genomic_position(int *flag, bwtint_t index, genomic_location *loc, char *cigar, uint32_t *mapq, bwtint_t seq_len,
                     int unmapped, int reverse, bwtint_t *offset, bwtint_t offInd) {
    if (index < 0 || index >= seq_len)
        *flag |= unmapped;
    if (((*flag) & unmapped) == 0) {
        if (index >= seq_len / 2) {
            index = (seq_len / 2) - (index - (seq_len / 2)) - cigar_len(cigar);
            *flag |= reverse;
        }
        int chr = binary_search(0, offInd - 1, offset, index);
        if (chr < 0) (*flag) |= unmapped;
        else {
            index -= offset[chr];
            loc->pos = index + 1;
            loc->rname = name[chr];
        }
    }
    if ((*flag) & unmapped) {
        mapq = 0;
        loc->pos = 0;
        loc->rname = star;
    }
}

int sam_generator(char *buffer, char *qname, int flag, uint32_t mapq, bwtint_t index, char *cigar, bwtint_t index2,
                  char *cigar2, ubyte_t *seq, ubyte_t *quality, int len, int len2, bwt_t *bwt, bwtint_t *offset,
                  bwtint_t offInd) {
    if (!cigar) cigar = star;
    if (!cigar2) cigar2 = star;
    genomic_location first, next;
    long long tlen = 0; // Template length: the size of fragment for paired-end reads
    char seq_string[len + 1], quality_string[len + 1];
    uint32_t tmp;
    get_genomic_position(&flag, index, &first, cigar, &mapq, bwt->seq_len, sf_unmapped, sf_reverse, offset, offInd);
    if (flag & sf_paired) {
        get_genomic_position(&flag, index2, &next, cigar2, &tmp, bwt->seq_len, sf_mate_unmapped, sf_mate_reverse,
                             offset, offInd);
        if ((flag & sf_unmapped) && ((flag & sf_mate_unmapped) == 0)) {
            // According to SAM format, when the mate is mapped but not read, the read rname and pos should be the same as mate
            first.rname = next.rname;
            first.pos = next.pos;
        } else if (((flag & sf_unmapped) == 0) && (flag & sf_mate_unmapped)) {
            next.rname = first.rname;
            next.pos = first.pos;
        }
    } else {
        next.pos = 0;
        next.rname = star;
    }
    if (flag & sf_pair_mapped) {
        if (first.pos < next.pos) tlen = next.pos + len2 - first.pos;
        else tlen = -(first.pos + len - next.pos);
    }
    if ((flag & sf_paired) && (flag & sf_unmapped) && ((flag & sf_mate_unmapped) == 0)) {
        // According to SAM format, when the mate is mapped but not read, the read rname and pos should be the same as mate
        first.rname = next.rname;
        first.pos = next.pos;
    }
    if ((flag & sf_paired) && (flag & sf_unmapped) && ((flag & sf_mate_unmapped) == 0)) {
        // According to SAM format, when the mate is mapped but not read, the read rname and pos should be the same as mate
        first.rname = next.rname;
        first.pos = next.pos;
    }
    if ((flag & sf_paired) && (strcmp(first.rname, next.rname) == 0) && strcmp(first.rname, "*")) next.rname = equal;

    seq_string[len] = quality_string[len] = 0;
    int i;
    if (flag & sf_not_primary) { // To make SAM file shorter, according to the SAM format requirements
        strcpy(seq_string, "*");
        strcpy(quality_string, "*");
    } else {
        if (!(flag & sf_reverse)) {
            for (i = 0; i < len; i++) {
                if (seq[i] < 4)
                    seq_string[i] = int_to_bp[seq[i]];
                else
                    seq_string[i] = 'N';
                quality_string[i] = quality[i];
            }
        } else {
            for (i = 0; i < len; i++) {
                if (seq[i] < 4)
                    seq_string[len - i - 1] = int_to_bp[3 - seq[i]];
                else
                    seq_string[len - i - 1] = 'N';
                quality_string[len - i - 1] = quality[i];
            }
            reverse_cigar(cigar);
        }
    }
    if (flag & sf_unmapped) cigar = star;
    return snprintf(buffer, MAX_SAM_LINE, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%"PRIu64"\t%lld\t%s\t%s\n", qname, flag,
                    first.rname, first.pos, mapq, cigar, next.rname, next.pos, tlen, seq_string, quality_string);
}

int sam_headers(char *buffer, bwtint_t *offset, int size, int buf_size) {
    int head = 0;
    int i;
    //fprintf(stderr, "offset[size] = %llu, offset[size - 1] = %llu, last = %llu, size = %d\n", offset[size], offset[size - 1], offset[size] - offset[size - 1], size);
    head += snprintf(buffer, buf_size - head, "@HD\tVN:1.0\tSO:unsorted\n");
    //fprintf(stderr, "buffer = %s offset[0] = %llu\n name[0] = %s\n", buffer, offset[0], name[0]);
    for (i = 0; i < size; i++) { //TODO
        head += snprintf(buffer + head, buf_size - head, "@SQ\tSN:%s\tLN:%"PRIu64"\n", name[i],
                         offset[i + 1] - offset[i]); //TODO
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


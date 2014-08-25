//
//  calcute_acuracy.c
//  
//
//  Created by Maryam Rabiee on 8/21/14.
//
//

#include <stdio.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    
    FILE *samFile;
    samFile = fopen("a.sam", "r");
    long long badAlignedReads = 0;
    long long notAlignedReads = 0;
    long long readNum = 0;
    char line[1000];
    int header = 1;
    char *qname, *rnext, *pnext, *seq_string, *quality_string,*rname, *cigar;
    int flag, i;
    uint64_t pos;
    uint32_t mapq;
    long long int tlen;
    rname = malloc(100 * sizeof(char));
    cigar = malloc(200 * sizeof(char));
    qname = malloc(100 * sizeof(char));
    rnext = malloc(100 * sizeof(char));
    pnext = malloc(100 * sizeof(char));
    seq_string = malloc(1000 * sizeof(char));
    quality_string = malloc(500 * sizeof(char));

    int stop = 0;
    while (1) {
            if (fgets(line, 1000, samFile) == NULL) 
                break;
            while (line[0] == '@'){
                if(fgets(line, 1000, samFile) == NULL){
                    stop = 1;
                    break;
                }
            }
            if(stop)
                break;
            readNum++;
            sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);
            char* tokens=strtok(qname, ":");
            tokens = strtok(NULL, ":");
    
            char *first,*second;
            first = strtok(tokens, "-");
            second = strtok(NULL, "-");
    
            printf( "f: %s   sec:%s \n",first,second  );
            uint64_t start, end;
            start = strtoll(first,NULL,10);
            end = strtoll(second,NULL,10);
            printf( "first: %" PRIu64 "   sec:%" PRIu64 " \n",start,end  );
            if(start-20 <= pos && pos <=end+20);
                //exact aligning
            else if(!strstr(cigar,"*"))
                badAlignedReads += 1;
    //readCigar(cigar, pos, seq_string, i,start ,end);

            fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n \n",qname, flag, rname, pos,mapq, cigar,rnext,pnext, tlen,seq_string,quality_string);

            fprintf(stdout,"\n cigar %s \n ",cigar);
            if(strstr(cigar,"*"))
                notAlignedReads++;
    }
    float accuracy = ((float)badAlignedReads/readNum)*100.0 ;
    fprintf(stdout, "number of not aligned reads: %lld \n number of reads with wrong alignment : %lld \n total reads : %lld \n percentage of bad aligned reads :%10f \n",notAlignedReads, badAlignedReads, readNum , accuracy);
    fclose(samFile);
    
}
void readCigar(char * cigar, uint64_t ref_i, char *seq_string, long readNum ) {
	//fprintf(stderr, "salam\n");
	int pos = 0;
	int value = 0;
	uint64_t ref_index = ref_i;
	long read_index = 0;
	char alignType;
	//printf("   %s\n", seq_string);
    //printf("cigar:   %s\n", cigar);
	while (1) {
		if (!isdigit(cigar[pos])) {
			//printf("   1salam\n");
			if (value > 0) {
                
				//printf("value:   %d",value);
				if (cigar[pos] == 'm') {
					int j;
                    
				} else if (cigar[pos] == 'd') {
					ref_index += value;
				} else if (cigar[pos] == 'i')
					read_index += value;
				else {
                    //					printf("*");
                    //readPenalties[readNum] += LONG_MAX;
					break;
				}
                
			}
            //			printf("*");
			if (cigar[pos] == 0)
				break;
			value = 0;
		} else {
			value = value * 10 + cigar[pos] - '0';
		}
		pos++;
	}
}


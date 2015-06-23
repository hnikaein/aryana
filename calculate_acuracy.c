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
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    printf("heyy\n");
    int each_wrong[4],i;
    for (i = 0 ; i<4; i++) {
        each_wrong[i]=0;
    }
    if ( argc != 2 ) {
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s filename", argv[0] );
    }
    else
    {
        // We assume argv[1] is a filename to open
        FILE *samFile = fopen( argv[1], "r" );
        /* fopen returns 0, the NULL pointer, on failure */
        if ( samFile == 0 )
            printf( "Could not open file\n" );
        printf("heyy");
        long count = 0;
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
        rname = (char*)malloc(100 * sizeof(char));
        cigar = (char*)malloc(200 * sizeof(char));
        qname = (char*)malloc(100 * sizeof(char));
        rnext = (char*)malloc(100 * sizeof(char));
        pnext = (char*)malloc(100 * sizeof(char));
        seq_string = (char*)malloc(1000 * sizeof(char));
        quality_string = (char*)malloc(500 * sizeof(char));

        int stop = 0;
        while (1) {
            if (fgets(line, 1000, samFile) == NULL)
                break;
            while (line[0] == '@') {
                if(fgets(line, 1000, samFile) == NULL) {
                    stop = 1;
                    break;
                }
            }
            if(stop)
                break;
            readNum++;
            printf("%d\n ",readNum);
            //1000000_chr21:30778456-30778555
            //968_>chr2:24442257-24442356_-o
            //20|chr3:75909157-75909256|+o

            char *copy = (char *)malloc(strlen(qname) + 1);
            strcpy(copy, qname);
            printf("heyy1\n");
            char* token=strtok(copy, ":");
            char* chrom = strtok(token, "|");
            chrom = strtok(NULL, "|");
            printf("hey2y\n");
            char *tokens = strtok(qname, ":");
            tokens = strtok(NULL, ":");
            printf( "chrom: %s   tokens:    %s \n",chrom,tokens  );
            printf("heyy3\n");
            char *copy2 = (char *)malloc(strlen(tokens) + 1);
            strcpy(copy2, tokens);
            printf("heyy4\n");
            char *first,*second,*type;
            first = strtok(tokens, "-");
            second = strtok(NULL, "-");
            second = strtok(second, "_");
            second = strtok(second, "|");

            type = strtok(copy2, "|");
            type = strtok(NULL, "|");

            fprintf(stdout, "%s   %s   %s\n",first,second,type);


            uint64_t start, end;
            start = strtoll(first,NULL,10);
            end = strtoll(second,NULL,10);
            //printf( "first: %" PRIu64 "   sec:%" PRIu64 " \n",start,end  );
            if(start-20 <= pos && pos <=end+20 && !strcmp(chrom,rname ));
            //exact aligning
            else if(!strstr(cigar,"*")) {
                badAlignedReads += 1;
                if (!strcmp(type, "+o"))
                    each_wrong[0]++;
                else if (!strcmp(type, "-o"))
                    each_wrong[1]++;
                else if (!strcmp(type, "+p"))
                    each_wrong[2]++;
                else if (!strcmp(type, "-p"))
                    each_wrong[3]++;

                //fprintf(stdout,"%s\n",line);
            }

            //fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n \n",qname, flag, rname, pos,mapq, cigar,rnext,pnext, tlen,seq_string,quality_string);

            if(strstr(cigar,"*"))
                notAlignedReads++;
            //      chromNum=0;
            //       alignedChrNum=0;
        }
        float accuracy = ((float)badAlignedReads/readNum)*100.0 ;
        float accuracy2 = ((float)notAlignedReads/readNum)*100.0 ;
        fprintf(stdout, "number of not aligned reads: %lld \n number of reads with wrong alignment : %lld \n total reads : %lld \n percentage of incorrectly aligned reads :%10f \n percentage of not aligned reads :%10f \n\n",notAlignedReads, badAlignedReads, readNum , accuracy,accuracy2);

        fprintf(stdout, "number of wrong alignment in each read type:\n +o : %d\n -o : %d\n +p : %d\n -p : %d \n",each_wrong[0],each_wrong[1],each_wrong[2],each_wrong[3]);
        fclose(samFile);
    }
}

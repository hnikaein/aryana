//
//  new_accuracy.cpp
//  aryana
//
//  Created by Maryam Rabiee on 9/30/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//

#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>


#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {


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

        //    FILE *samFile;
        //    samFile = fopen(samName, "r");
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

        rname = new char[100];
        cigar = new char[200];
        qname = new char[100];
        rnext = new char[100];
        pnext = new char[100];
        seq_string = new char[1000];
        quality_string = new char[500];

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

            sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos,&mapq, cigar,rnext,pnext, &tlen,seq_string,quality_string);

            string readname(seq_string);
            //1000000_chr21:30778456-30778555
            //968_>chr2:24442257-24442356_-o

            char *copy = (char *)malloc(strlen(qname) + 1);
            strcpy(copy, qname);

            char* token=strtok(copy, ":");
            char* chrom = strtok(token, ">");
            chrom = strtok(NULL, ">");

            char *tokens = strtok(qname, ":");
            tokens = strtok(NULL, ":");
            //printf( "chrom: %s   tokens:    %s \n",chrom,tokens  );

            char *first,*second;
            first = strtok(tokens, "-");
            second = strtok(NULL, "-");

            second = strtok(second, "_");

            uint64_t start, end;
            start = strtoll(first,NULL,10);
            end = strtoll(second,NULL,10);
            //printf( "first: %" PRIu64 "   sec:%" PRIu64 " \n",start,end  );
            if(start-20 <= pos && pos <=end+20 && !strcmp(chrom,rname ));
            //exact aligning
            else if(!strstr(cigar,"*")) {
                badAlignedReads += 1;

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
        fprintf(stdout, "number of not aligned reads: %lld \n number of reads with wrong alignment : %lld \n total reads : %lld \n percentage of bad aligned reads :%10f \n percentage of not aligned reads :%10f \n",notAlignedReads, badAlignedReads, readNum , accuracy,accuracy2);
        fclose(samFile);
    }
}











































void computeMethylation() {

    if(lines.size() == 0)
        return;
    int lastchecked = lines[0].pos;
    if(cytosines.size()!=0) {
        lastchecked= cytosines[cytosines.size()-1].pos+1;
        if (lines[0].pos > lastchecked || strcmp(cytosines[cytosines.size()-1].chr,lines[0].chr)) {//if chrom has changed or start of the read has passed lastchecked pos
            lastchecked = lines[0].pos;
        }
    }
    long refPos = chrom[ChromIndex(lines[0].chr)].chrStart;

    for(long i = lastchecked ; i< (lines[0].seq_string.size()+lines[0].pos ) ; i++) {
        if(toupper(reference[refPos + i-1]) == 'C') {
            setPointer(i , lines[0].chr);
            cytosines.push_back(Cytosine(i ,lines[0].chr,lines[0].strand));
            for(int j=0 ; j< secondPointer ; j++) {
                if( lines[j].strand == '+') {
                    int relative_pos = lines[j].pos + lines[j].seq_string.size();
                    if(i < relative_pos) {
                        if (lines[j].seq_string[i - lines[j].pos] == 'C') {
                            cytosines[cytosines.size()-1].methylated++;
                        }
                        else if(lines[j].seq_string[i - lines[j].pos] == 'T')
                            cytosines[cytosines.size()-1].unmethylated++;
                    }
                }
            }
        }
        else if(toupper(reference[refPos + i-1]) == 'G') {
            setPointer(i , lines[0].chr);
            cytosines.push_back(Cytosine(i ,lines[0].chr,lines[0].strand));
            for(int j=0 ; j< secondPointer ; j++) {
                if(lines[j].strand == '-') {

                    if(!(lines[j].pos < (lines[0].pos+lines[0].seq_string.size())))
                        break;
                    int relative_pos = lines[j].pos + lines[j].seq_string.size();
                    if(i < relative_pos ) {

                        if (lines[j].seq_string[i - lines[j].pos] == 'G')
                            cytosines[cytosines.size()-1].methylated++;

                        else if(lines[j].seq_string[i - lines[j].pos] == 'A')
                            cytosines[cytosines.size()-1].unmethylated++;
                    }
                }
            }
        }
    }

    lines.erase(lines.begin());
    secondPointer--;
    computeMethylation();

}
void setPointer(int pos, char * chr) {
    while(secondPointer < lines.size() && lines[secondPointer].pos <= pos && !strcmp(lines[secondPointer].chr,chr)) {
        secondPointer++;
    }
    if (secondPointer == lines.size() && lines[secondPointer].pos <= pos && !strcmp(lines[secondPointer].chr,chr)) {
        //int temp = i;
        int res = readSamFile(samFile);
        if(res != -1)
            setPointer(pos, chr);

    }
    else if(secondPointer == lines.size())
        secondPointer--;
}

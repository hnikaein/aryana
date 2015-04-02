//
//  methylation_state.c
//  aryana
//
//  Created by Maryam Rabiee on 9/15/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//

#include <stdio.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>


#define maxGenomeSize 1e5
#define maxChromosomeNum 1000
#define READ_LENGHT 100
unsigned long gs;
char chromName[maxChromosomeNum][100];
int chromNum;

char * reference;
uint32_t reference_size;
uint32_t reference_reminder;

struct line {
    char line[1000];
    char chr[30];
    uint64_t start, end , pos;
    char *qname, *seq_string,*rname, *cigar;

};
struct line lines[1000];

struct Chrom
{
    char * chrName;
    long chrStart;
    long chrNum;
};

struct Chrom chrom[maxChromosomeNum] ;

int main(int argc, char *argv[]) {
    char *referenceName, *annotationFile;

    if ( argc != 2 ) {
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s filename", argv[0] );
    }
    else
    {
        ReadGenome(referenceName);
        // We assume argv[1] is a filename to open
        FILE *samFile = fopen( argv[1], "r" );
        /* fopen returns 0, the NULL pointer, on failure */
        if ( samFile == 0 )
            printf( "Could not open file\n" );

        long long readNum = 0;
        char line[1000];
        int header = 1;
        char *qname, *rnext, *pnext, *seq_string, *quality_string,*rname, *cigar;
        int flag, i;
        uint64_t pos;
        uint32_t mapq;
        long long int tlen;
        for(i=0 ; i<1000; i++) {
            lines[i].cigar = malloc(200 * sizeof(char));
            lines[i].rname = malloc(100 * sizeof(char));
            lines[i].seq_string = malloc(1000 * sizeof(char));
        }
        rname = malloc(100 * sizeof(char));
        cigar = malloc(200 * sizeof(char));
        qname = malloc(100 * sizeof(char));
        rnext = malloc(100 * sizeof(char));
        pnext = malloc(100 * sizeof(char));
        seq_string = malloc(1000 * sizeof(char));
        quality_string = malloc(500 * sizeof(char));

        int count_to_tousand=0;
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

            sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, lines[count_to_tousand].chr, &lines[count_to_tousand].pos ,&mapq, lines[count_to_tousand].cigar,rnext,pnext, &tlen,lines[count_to_tousand].seq_string,quality_string);
            count_to_tousand++;
            //1000000_chr21:30778456-30778555
            //968_>chr2:24442257-24442356

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

            //uint64_t start, end;
            lines[count_to_tousand].start = strtoll(first,NULL,10);
            lines[count_to_tousand].end = strtoll(second,NULL,10);
            if (count_to_tousand==1000) {
                computeMethylation(lines[count_to_tousand].start,lines[count_to_tousand].end , lines[count_to_tousand].seq_string);
                count_to_tousand=0;
            }


        }//while
        fclose(samFile);
    }

}//main


int ReadGenome(char * genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info;
    if (stat(genomeFile, &file_info) == -1) {
        fprintf(stderr,
                "Could not get the information of file %s\nplease make sure the file exists\n",
                genomeFile);
        return -1;
    }
    FILE *fp;
    fp = fopen(genomeFile, "r");
    off_t file_size_bytes = file_info.st_size;
    reference_size = ceil(((double) file_size_bytes) / (double) (sizeof(char)));
    reference = (char *) malloc(reference_size * sizeof(char));
    memset(reference, 0, reference_size * sizeof(char));
    gs = 0;
    chromNum = 0;
    //fprintf(stderr, "Reading genome...\n");
    char fLine[10000];
    while (! feof(fp)) {
        // fprintf(stderr, "Reading genome1...\n");
        int n = fscanf(fp, "%s\n", fLine);
        if (n == EOF) break;
        n = strlen(fLine);
        // fprintf(stderr, "Reading genome2 :n %d ...\n",n);
        if (fLine[0] == '>') {
            //chromPos[chromNum++] = gs;
            chrom[chromNum++].chrStart = gs;
            char * temp = fLine;
            temp++;
            int lenght = strlen(temp);
            chrom[chromNum-1].chrName = (char *) malloc(lenght * sizeof(char));
            strcpy(chrom[chromNum-1].chrName, temp);
            //fprintf(stderr, "Reading genome3...\n");
            if (chromNum >= 1)
                fprintf(stderr, "chrName: %s, chrStart: %ld\n",chrom[chromNum-1].chrName, chrom[chromNum-1].chrStart);
            fprintf(stderr, "%s\n",fLine);
            //strcpy(chromName[chromNum - 1], fLine);
        } else {
            //            ToUpper(fLine);
            memcpy(reference+gs, fLine, n);
            gs += n;
        }
    }
    //chromPos[chromNum] = gs;
    // fprintf(stderr, "Lenght: %ld\n",chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(fp);
}


long methylated = 0;
void computeMethylation() {
    int j=0;
    int lastchecked ;

    while(j < 1000) {


        int pos = lines[j].pos;
        char *chrom = lines[j].chr;
        pch=strchr(lines[j].seq_string,'C');
        int refPos2 = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart;
        int refPos = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart +
                     pch - lines[j].seq_string + 1;
        lastchecked = j;
        if(reference[refPos] == 'C' && reference[refPos+1] == 'G') {
            methylated ++;
        }
        else {
            outOfcontext++;
            pch=strchr(lines[j].seq_string,'C');
            if(pch == NULL) {
                j = lastchecked + 1;
                continue;
            }

        }
        while (!strcmp(line[j].chr, chrom) && (pos + READ_LENGHT) >= lines[j].pos && j < 1000) {
            char * pch;
            pch2=strchr(lines[j].seq_string,'C');               ////////////handle lowercase too
            if(pch2 =pch)
                methylfolan++;
            else
                j++;

//            int refPos2 = lines[j].pos + chrom[ChromIndex(lines[j].rname)].chrStart;
//            int loc = pch - lines[j].seq_string+1;


        }

    }

}

void computeMethylation2(char *chr) {
    int i;

    int index = ChromIndex(chr);
    int start = chrom[index].chrStart;
    int end = chrom[index].
}
int ChromIndex(char * chr) {
    int i = 0;
    for(i; i < maxChromosomeNum; i++) {
        if(chrom[i].chrName == NULL)
            break;
        if(strcmp(chrom[i].chrName, chr) == 0)
            return i;
    }
    return -1;

}

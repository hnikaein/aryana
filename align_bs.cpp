#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <libgen.h>
#include "utils.h"
#include "bwt.h"
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <iostream>
using namespace std;

#define maxChromosomeNum 1000
#define numberOfGenomes 5
#define maxReadLength 2000
#define maxChrNameLength 100
#define maxSamFileLineLength 10000
//const long maxGenomeSize = 4e9;
//const int maxChromosomeNum = 1000;

long long chromPos[maxChromosomeNum];
unsigned long long gs;
char chromName[maxChromosomeNum][maxChrNameLength];
char * genome;
int chromNum;
char * genomeFile, *annotationFile, *outputFile, *samFilePath, *referenceName;

struct Island {
    int chr;
    long long start, end;
    Island (int c, long long s, long long e) {
        chr = c;
        start = min(s,e);
        end = max(s,e);
    }
    bool operator < (const Island & val) const {
        return chr < val.chr || (chr == val.chr && end < val.start);
    }
};

vector <Island> islands;
long islandsNum = 0;
int highPenalty = 0, medPenalty = 0, lowPenalty = 0;
long penalties;
long long readNum;
long readPenalties[numberOfGenomes];

char *samNames[numberOfGenomes];
FILE *samFiles[numberOfGenomes];
char str[maxReadLength];


struct Chrom
{
    char * chrName;
    long chrStart;
}; // for converting refrence position to chromosome based position

struct Chrom chrom[maxChromosomeNum] ;
char * reference;
uint32_t reference_size;
uint32_t reference_reminder;


// Reads genome and saves start position of chromosomes
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
    char fLine[10000];
    while (! feof(fp)) {
        int n = fscanf(fp, "%s\n", fLine);
        if (n == EOF) break;
        n = strlen(fLine);
        if (fLine[0] == '>') {
            chrom[chromNum++].chrStart = gs;
            char * temp = fLine;
            temp++;
            int lenght = strlen(temp);
            chrom[chromNum-1].chrName = (char *) malloc(lenght * sizeof(char));
            strcpy(chrom[chromNum-1].chrName, temp);
            if (chromNum >= 1)
                fprintf(stderr, "chrName: %s, chrStart: %ld\n",chrom[chromNum-1].chrName, chrom[chromNum-1].chrStart);
            fprintf(stderr, "%s\n",fLine);
        } else {
            memcpy(reference+gs, fLine, n);
            gs += n;
        }
    }
    chromPos[chromNum] = gs;
    fprintf(stderr, "Lenght: %lld\n",chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(fp);
	return 0;
}


void reverseRead(char *s) {
    char  c;
    char* s0 = s - 1;
    char* s1 = s;
    while (*s1) ++s1;
    while (s1-- > ++s0)
    {
        c   = *s0;
        *s0 = *s1;
        *s1 =  c;
    }
    
    return;
}

void complementRead(char *read) {
    int i;
    for (i=0; i< strlen(read); i++) {
        switch (read[i]) {
            case 'T':
            case 't':
                read[i]='A';
                break;
            case 'C':
            case 'c':
                read[i]='G';
                break;
            case 'G':
            case 'g':
                read[i]='C';
                break;
            case 'A':
            case 'a':
                read[i]='T';
                break;
            default:
                break;
        }
    }
    return;
}

// returns minimum penalty from scores for each genome , if penalties are equal , returns original genome
int min_penalty() {
    int j = 0;
    long min = readPenalties[0];
    for (int i = 1; i< numberOfGenomes ; i++){
       // cerr<<readPenalties[i]<<endl;
        if(min > readPenalties[i]) {
            min = readPenalties[i];
            j = i;
        }
    }
    return j;
}

char getNuc(uint64_t place) {
    return reference[place];
}

inline void ToLower(char * s) {
    int i = 0;
    for (i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z')
            s[i] += 'a' - 'A';
}


int ChromIndex(char * chr) {
    for(int i=0; i < maxChromosomeNum; i++) {
        if(chrom[i].chrName == NULL)
            break;
        if(strcmp(chrom[i].chrName, chr) == 0)
            return i;
    }
    return -1;
    
}


void ReadCpGIslands(char * annotationFile) {
    FILE* file = fopen(annotationFile, "r"); /* should check the result */
    char fLine[10000], chrom[maxChrNameLength];
    uint64_t start, end;
    do {
        if(NULL==fgets(fLine, sizeof(fLine), file)) {
            fclose(file);
            break;
        }
        sscanf(fLine, "%s %lld %lld", chrom, &start, &end);
        islands.push_back(Island(ChromIndex(chrom), start, end));
    } while (true);
    sort(islands.begin(), islands.end());
}

void setPenalties(int p1, int p2, int p3) {
    highPenalty = p1;
    medPenalty = p2;
    lowPenalty = p3;
    return;
}

bool isInIsland(int chr, long long pos) {
    return binary_search(islands.begin(), islands.end(), Island(chr, pos, pos));
}


void CalcPenalties(uint64_t ref_i, char read, long readNum, int chr,uint64_t chrPos,int flag,int flag2) {
    read = toupper(read);
    reference[ref_i] = toupper(reference[ref_i]);
    reference[ref_i+1] = toupper(reference[ref_i+1]);
    reference[ref_i-1] = toupper(reference[ref_i-1]);
    
    if(flag2) {
        if (read == 'A' || read == 'G'){
            if (getNuc(ref_i) != read)
                readPenalties[readNum] += highPenalty;
        }
            else { //read = C or T
                
                if (read == 'T' && getNuc(ref_i) == 'C') {
                    if (getNuc(ref_i+1) == 'G') { // in the CpG context
                        if (isInIsland(chr,chrPos)) // in CpG and also island
                            readPenalties[readNum] += 0;
                        else
                            readPenalties[readNum] += highPenalty;
                        
                    } else // out of CpG context
                        
                        readPenalties[readNum] += lowPenalty;
                    
                } else if (read == 'C' && getNuc(ref_i) == 'C') {
                    if (getNuc(ref_i+1) == 'G') { // in the CpG context
                        int temp = isInIsland(chr , chrPos);
                        if (temp == 1)  // in CpG and also island
                            readPenalties[readNum] += highPenalty;
                        else            // in CpG and out of island
                            readPenalties[readNum] += lowPenalty;
                        
                    } else              // out of CpG context
                        readPenalties[readNum] += highPenalty;
                    
                } else if (read == 'C' && getNuc(ref_i) == 'T')
                    readPenalties[readNum] += highPenalty;
                
            }
    }
    else {
        if (read == 'C' || read == 'T'){
            if (getNuc(ref_i) != read)
                readPenalties[readNum] += highPenalty;
        }
            else { //read = G or A
                
                if (read == 'A' && getNuc(ref_i) == 'G') {
                    if (getNuc(ref_i-1) == 'C') { // in the CpG context
                        if (isInIsland(chr , chrPos)) { // in CpG and also island
                            readPenalties[readNum] += 0;
                        } else {
                            readPenalties[readNum] += highPenalty;
                        }
                    } else { // out of CpG context
                        
                        readPenalties[readNum] += lowPenalty;
                    }
                } else if (read == 'G' && getNuc(ref_i) == 'G') {
                    if (getNuc(ref_i-1) == 'C') { // in the CpG context
                        int temp = isInIsland(chr , chrPos);
                        if (temp == 1) { // in CpG and also island
                            readPenalties[readNum] += highPenalty;
                        } else {
                            readPenalties[readNum] += lowPenalty;
                        }
                    } else { // out of CpG context
                        readPenalties[readNum] += highPenalty;
                    }
                    
                } else if (read == 'G' && getNuc(ref_i) == 'A')
                    readPenalties[readNum] += highPenalty;
            }
        
    }
    
}

// Reads Cigar sequence and call CalcPenalties for calculating penalties for each base
void readCigar(char * cigar, uint64_t ref_i, char *seq_string, long long readNum, int chr,uint64_t chrPos,int flag, int refNum) {
    int pos = 0;
    int value = 0;
    uint64_t ref_index = ref_i;
    long read_index = 0;
    strcpy(str, seq_string);
    while (1) {
        if (!isdigit(cigar[pos])) {
            if (value > 0) {
                if (cigar[pos] == 'm') {
                    int j;
                    for (j = 0; j < value; j++) {
                        CalcPenalties(ref_index, seq_string[read_index], readNum,chr,chrPos,flag,refNum);
                        ref_index++;
                        read_index++;
                    }
                } else if (cigar[pos] == 'd') {
                    ref_index += value;
                    readPenalties[readNum] += highPenalty * value;      //high penalty for insertion
                } else if (cigar[pos] == 'i') {
                    read_index += value;
                    readPenalties[readNum] += highPenalty * value;      //high penalty for deletion
                }
            }
            else if(cigar[pos]=='*') {
                readPenalties[readNum] += LONG_MAX;     // maximum penalty for not aligned reads
                break;
            }
            if (cigar[pos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[pos] - '0';
        }
        pos++;
    }
}



int main(int argc, char *argv[]) {
    string header="";
    char line[maxSamFileLineLength];
    if (argc < 6) {
        fprintf(stderr, "Need more inputs\n");
        return -1;
    }
    static struct option long_options[] = {
        { "samfiles", required_argument, 0, 's' },
        { "reference", required_argument, 0, 'x' },
        { "CpG-islands", required_argument, 0, 'c' },
        { "penalties", required_argument, 0, 'p' }
    };
    int option_index = 0;
    char* samFileName = 0, *tmp;
    int c, i;
    while ((c = getopt_long(argc, argv, "x:s:c:p:", long_options, &option_index))
           >= 0) {
        switch (c) {
            case 0:
                break;
            case 1:
                break;
            case 'x':
                referenceName = (char *) malloc(strlen(optarg)+1);
                strcpy(referenceName, optarg);
                break;
            case 's':
                samFileName = strdup(optarg);
                tmp = strdup(samFileName);
                samFilePath = dirname(tmp);
                for(i=0; i<numberOfGenomes; i++) {
                    samNames[i] = (char *) malloc(strlen(samFileName) + 10);
                    sprintf(samNames[i], "%s-%d", samFileName, i);
                }
                break;
            case 'c':
                annotationFile = (char *) malloc(strlen(optarg)+1);
                strcpy(annotationFile, optarg);
                break;
            case 'p':
                lowPenalty = atoi(optarg);
                medPenalty = atoi(argv[optind]);
                highPenalty = atoi(argv[optind + 1]);
                optind = optind + 2;
                break;
        }
    }
    if (! samFileName) {
        cerr << "Error: input SAM file name not found." << endl;
        return 1;
    }
    ReadGenome(referenceName);
    
    ReadCpGIslands(annotationFile);
    
    // Printing header of original samfile to standard output
    FILE *samFile = fopen(samNames[0], "r");
    fgets(line, 10000, samFile);
    while (line[0] == '@'){
        header += line;
        fgets(line, 10000, samFile);
    }
    fprintf(stderr,"%s",header.c_str());
    fclose(samFile);
    
    
    char command[] = "LC_ALL=C sort -k 1 ";
    char buf[strlen(samNames[0]) + 100];
    for (i = 0; i < numberOfGenomes; i++) {
        sprintf(buf, "%s%s > %s/samFile%d.sam", command, samNames[i], samFilePath, i+1);
        system(buf);
        sprintf(buf, "%s/samFile%d.sam", samFilePath, i+1);
        samFiles[i] = fopen(buf, "r");

    }
    
    char *rname[numberOfGenomes], *cigar[numberOfGenomes];
    int flag[numberOfGenomes];
    uint64_t pos[numberOfGenomes];
    uint32_t mapq[numberOfGenomes];
    long long int tlen;
    for (i = 0; i < numberOfGenomes; i++) {
        rname[i] = new char[maxChrNameLength];
        cigar[i] = new char[maxReadLength * 2];
    }
    char qname[100], rnext[100], pnext[100], seq_string[maxReadLength], quality_string[maxReadLength], copy[maxReadLength]; // (Ali) please double check the limitations for qname, rnext, pnext and copy
    int chosen[numberOfGenomes];
    int j=0;
    for (; j < numberOfGenomes; j++)
        chosen[j] = 0;
    
    int stop = 0;
    
    while (1 && !stop) {
        //int reversed = 0;///////////////////////////bugesh dorost she
        for (i = 0; i < numberOfGenomes && !stop; i++) {
            do{
                if(fgets(line, 10000, samFiles[i]) == NULL){
                    stop = 1;
                    break;
                }
            } while (line[0] == '@');
            if(stop)
                break;
            //reversed = 0;
            sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag[i], rname[i], &pos[i],&mapq[i], cigar[i],rnext,pnext, &tlen,seq_string,quality_string);
            int index = ChromIndex(rname[i]);
            if(index == -1){
                readPenalties[i] += LONG_MAX;
                continue;
            }
            int flag2 = 0;
            if(i <=2) // if read has been aligned to first four genomes
                flag2 = 1;
            if(i > 2) // if read has been aligned to last three genome
                flag2 = 0;
            if(flag[i] == 16){
                //reversed = 1;
                strcpy(copy,seq_string); // copy read for printing it later
                reverseRead(copy);
                complementRead(copy);
                readCigar(cigar[i], pos[i]+chrom[index].chrStart-1, copy, i ,index,pos[i],flag[i],flag2);

            }
            else
                readCigar(cigar[i], pos[i]+chrom[index].chrStart-1, seq_string, i ,index,pos[i],flag[i],flag2);
        }
        if(stop)
            break;

        int min = min_penalty();
        chosen[min]++; // shows how many times a genome has been selected
        fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\t%d\n",qname, flag[min], rname[min], pos[min],mapq[min], cigar[min],rnext,pnext, tlen,seq_string,quality_string, min);
        int j=0;
        for (; j < numberOfGenomes; j++)
            readPenalties[j]=0;
    }
    
    for (j=0; j < numberOfGenomes; j++)
        fprintf(stderr, "Number of reads aligned to genome %d: %d\n", j, chosen[j]);
    
}


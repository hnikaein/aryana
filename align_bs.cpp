#define __STDC_FORMAT_MACROS
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <cmath>
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

#define numberOfGenomes 5
#define maxReadLength 10000
#define maxChrNameLength 10000
#define maxSamFileLineLength 20000
const int unalignedPenalty = 1000000;
const int discordPenalty = 100000;
unsigned long long gs;
vector <string> chromName;
char * genome;
int chromNum, minD = 0, maxD = 10000, limit = 5;
char * genomeFile, *annotationFile, *samFilePath, *referenceName;
FILE * outputFile = NULL, * headerFile = NULL;
bool paired = false;
int chosen[numberOfGenomes];
char line[numberOfGenomes][maxSamFileLineLength], line2[numberOfGenomes][maxSamFileLineLength];

enum orientation_t {fr, rf, ff, all} orientation = all;

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
    string chrName;
    long long chrStart;
}; // for converting refrence position to chromosome based position

vector <Chrom> chrom;


// Reads genome and saves start position of chromosomes
int ReadGenome(char * genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info;
    FILE *f;
    if (stat(genomeFile, &file_info) == -1) {
        fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", genomeFile);
        return -1;
    }

    f = fopen(genomeFile, "r");
    if (! f) {
        fprintf(stderr, "Error opening reference file: %s\n", genomeFile);
        exit(-1);
    }
    off_t file_size_bytes = file_info.st_size;
    long long reference_size = 1 + ceil(((double) file_size_bytes) / (double) (sizeof(char)));
    genome = (char *) malloc(reference_size * sizeof(char));
    memset(genome, 0, reference_size * sizeof(char));
    gs = 0;
    chromNum = 0;
    fprintf(stderr, "Reading genome...\n");
    char fLineMain[1000000];
    Chrom ch;
    while (! feof(f)) {
        if (! fgets(fLineMain, sizeof(fLineMain), f)) break;
        int n = strlen(fLineMain), start = 0;
        while (n > 0 && fLineMain[n-1] <= ' ') n--;
        fLineMain[n] = 0;
        while (start < n && fLineMain[start] <= ' ') start++;
        if (start >= n) continue;
        char * fLine = fLineMain + start;
        n -= start;

        if (fLine[0] == '>') {
            ch.chrStart = gs;
            string name = fLine;
            if (name.find(" ") != string::npos) name = name.substr(1, name.find(" ")-1);
            else name = name.substr(1, name.size() - 1);
            //cerr << name << endl;
            ch.chrName = name;
            chrom.push_back(ch);
            chromNum++;
        } else {
            memcpy(genome+gs, fLine, n);
            gs += n;
        }
    }
    fclose(f);
    return 0;
}


// returns minimum penalty from scores for each genome , if penalties are equal , returns original genome
int min_penalty() {
    int j = 0;
    long min = readPenalties[0];
    for (int i = 1; i< numberOfGenomes ; i++) {
        // cerr<<readPenalties[i]<<endl;
        if(min > readPenalties[i]) {
            min = readPenalties[i];
            j = i;
        }
    }
    return j;
}

char getNuc(uint64_t place) {
    return genome[place];
}

inline void ToLower(char * s) {
    unsigned int i = 0;
    for (i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z')
            s[i] += 'a' - 'A';
}


int ChromIndex(char * chr) {
    for(unsigned int i=0; i < chrom.size(); i++) {
        if(chrom[i].chrName == chr)
            return i;
    }
    return -1;

}


void ReadCpGIslands(char * annotationFile) {
    FILE* file = fopen(annotationFile, "r"); /* should check the result */
    char fLine[10000], chrom[maxChrNameLength];
    unsigned long long start, end;
    do {
        if(NULL==fgets(fLine, sizeof(fLine), file)) {
            fclose(file);
            break;
        }
        sscanf(fLine, "%s %llu %llu", chrom, &start, &end);
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
    genome[ref_i] = toupper(genome[ref_i]);
    genome[ref_i+1] = toupper(genome[ref_i+1]);
    genome[ref_i-1] = toupper(genome[ref_i-1]);

    if(flag2) {
        if (read == 'A' || read == 'G') {
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
        if (read == 'C' || read == 'T') {
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
                readPenalties[readNum] += unalignedPenalty;     // maximum penalty for not aligned reads
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

void WriteHeader() {
    // Printing header of original samfile to standard output
    char line[maxSamFileLineLength];
    string header;
    FILE *samFile = fopen(samNames[0], "r");
    fgets(line, 10000, samFile);
    while (line[0] == '@') {
        header += line;
        fgets(line, 10000, samFile);
    }
    fprintf(headerFile,"%s",header.c_str());
    fclose(samFile);
    if (headerFile != outputFile) fclose(headerFile);
}

bool ReadLine(char * line, int len, FILE * f) {
    do {
        if(fgets(line, len, f) == NULL) return false;
    } while (line[0] == '@');
    return true;
}

void Process() {
    for (int i = 0; i < numberOfGenomes; i++) samFiles[i] = fopen(samNames[i], "r");
    char *rname[numberOfGenomes], *cigar[numberOfGenomes],*rname2[numberOfGenomes], *cigar2[numberOfGenomes];
    int flag[numberOfGenomes],flagTwo[numberOfGenomes];
    uint64_t pos[numberOfGenomes],pos2[numberOfGenomes];
    uint32_t mapq[numberOfGenomes],mapq2[numberOfGenomes];
    long long int tlen,tlen2;
    for (int i = 0; i < numberOfGenomes; i++) {
        rname[i] = new char[maxChrNameLength];
        cigar[i] = new char[maxReadLength * 2];
        rname2[i] = new char[maxChrNameLength];
        cigar2[i] = new char[maxReadLength * 2];
    }
    char qname[10000], rnext[10000], pnext[10000], seq_string[maxReadLength], quality_string[maxReadLength]; // (Ali) please double check the limitations for qname, rnext, pnext and copy
    char qname2[10000], rnext2[10000], pnext2[10000], seq_string2[maxReadLength], quality_string2[maxReadLength];
    bzero(chosen, sizeof(chosen));

    while (true) {
        for (int i = 0; i < numberOfGenomes; i++) {
            readPenalties[i] = 0;
            if (! ReadLine(line[i], maxSamFileLineLength, samFiles[i])) return;
            sscanf(line[i], "%s\t%d\t%s\t%" PRIu64 "\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n" ,qname, &flag[i], rname[i], &pos[i],&mapq[i], cigar[i],rnext,pnext, &tlen,seq_string,quality_string);
            int index = ChromIndex(rname[i]);
            if(index == -1) readPenalties[i] += unalignedPenalty;
            else {
                int flag2 = 0;
                if(i <=2) // if read has been aligned to first four genomes
                    flag2 = 1;
                if(i > 2) // if read has been aligned to last three genome
                    flag2 = 0;
                readCigar(cigar[i], pos[i]+chrom[index].chrStart-1, seq_string, i ,index,pos[i],flag[i],flag2);
            }

            if (paired) {
                if (! ReadLine(line2[i], maxSamFileLineLength, samFiles[i])) {
                    cerr << "Error: unexpected end of input SAM file while reading the paired end.\n";
                    return;
                }
                sscanf(line2[i],"%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname2, &flagTwo[i], rname2[i], &pos2[i],&mapq2[i], cigar2[i],rnext2,pnext2, &tlen2,seq_string2,quality_string2);
                int index2 = ChromIndex(rname2[i]);
                if(index2 == -1) readPenalties[i] += unalignedPenalty;
                else {
                    int flag2 = 0;
                    if(i <=2) // if read has been aligned to first four genomes
                        flag2 = 1;
                    if(i > 2) // if read has been aligned to last three genome
                        flag2 = 0;
                    readCigar(cigar2[i], pos2[i]+chrom[index2].chrStart-1, seq_string2, i ,index2,pos2[i],flagTwo[i],flag2);
                }
                orientation_t orien;
                if (index >= 0 && index2 >= 0) {
                    if ((flag[i] & 16) == (flagTwo[i] & 16)) orien = ff;
                    else if (((flag[i] & 16) == 0 && pos[i] < pos2[i]) || ((flag[i] & 16) != 0 && pos[i] > pos2[i])) orien = fr;
                    else orien = rf;
                    if (index != index2 || llabs((signed) (pos[i] - pos2[i])) < minD || llabs((signed)(pos[i] - pos2[i])) > maxD || (orientation != all && orientation != orien)) readPenalties[i] += discordPenalty;
                }
            }
        }
        int min = min_penalty();
        chosen[min]++; // shows how many times a genome has been selected
        fprintf(outputFile, "%s", line[min]);
        if (paired) fprintf(outputFile, "%s", line2[min]);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr, "Usage: align_bs -x <reference genome> -s <input sam files> -c <CpG-island map> -p <penalties> [-o <output SAM file>] [-h <SAM header file>] [-l/--limit <maximum allowed mismatches>]\n");
        fprintf(stderr, "Additional arguments for paired-end reads:\n");
        fprintf(stderr, "    -P to consider the input as paired-end, -m/--min for minimum and -M/--max for maximum distance, and --fr/rf/ff for orientation.\n");
        return -1;
    }
    static struct option long_options[] = {
        { "samfiles", required_argument, 0, 's' },
        { "reference", required_argument, 0, 'x' },
        { "CpG-islands", required_argument, 0, 'c' },
        { "penalties", required_argument, 0, 'p' },
        { "output", required_argument, 0, 'o' },
        { "header", required_argument, 0, 'h' },
        { "paired-end", no_argument, 0, 'U' },
        { "min", required_argument, 0, 'm'},
        { "max", required_argument, 0, 'M'},
        { "fr", no_argument, 0, 1},
        { "rf", no_argument, 0, 2},
        { "ff", no_argument, 0, 3},
        { "limit", required_argument, 0, 'l'}
    };
    int option_index = 0;
    char* samFileName = 0, *tmp, c;
    outputFile = stdout;
    while ((c = getopt_long(argc, argv, "x:s:c:p:o:h:Pm:M:123l:", long_options, &option_index)) >= 0) {
        switch (c) {
        case 'x':
            referenceName = strdup(optarg);
            break;
        case 's':
            samFileName = strdup(optarg);
            tmp = strdup(samFileName);
            samFilePath = dirname(tmp);
            for(int i=0; i<numberOfGenomes; i++) {
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
        case 'o':
            outputFile = fopen(optarg, "w");
            break;
        case 'h':
            headerFile = fopen(optarg, "w");
            break;
        case 'P':
            paired = true;
            break;
        case 'm':
            minD = atoi(optarg);
            break;
        case 'M':
            maxD = atoi(optarg);
            break;
        case 1:
            orientation = fr;
            break;
        case 2:
            orientation = rf;
            break;
        case 3:
            orientation = ff;
            break;
        case 'l':
            limit = atoi(optarg);
            break;
        }
    }
    if (! headerFile) headerFile = outputFile;

    if (! samFileName) {
        cerr << "Error: input SAM file name not found." << endl;
        return 1;
    }
    ReadGenome(referenceName);
    ReadCpGIslands(annotationFile);
    WriteHeader();
    Process();
    for (int j=0; j < numberOfGenomes; j++)
        fprintf(stderr, "Number of reads aligned to genome %d: %d\n", j, chosen[j]);

    if (outputFile != stdout) fclose(outputFile);
}


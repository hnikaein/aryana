//
//  methylation.cpp
//  aryana
//
//  Created by Maryam Rabiee on 9/20/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//


#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


#define READ_SAM_BY 10
#define maxReadLength 2000
#define BUFFER_SIZE (maxReadLength + 1) // Should be higher than maximum read length
#define maxChrNameLength 50
#define maxReadNameLength 100
#define maxSamFileLineLength 10000
unsigned long long gs;
long long amb = 0; // the number of ambiguous reads, for which CT conversions equals GA conversions

int chromNum;
int debug=0;

char * genome;
FILE *samFile;
bool allCytosines = false; // Whether or not print the information of cytosines neither methylated nor in CpG context

void log(char* s) {
#ifdef DEBUG
    cout << s << endl;
#endif
}

struct SamRecord {
    char line[maxReadLength];
    char chr[maxChrNameLength];
    uint64_t start , pos;
    char cigar[maxReadLength * 2], cigar2[maxReadLength * 2];
    char rname[maxReadNameLength];
    string seq_string;
    char strand;
} line; // (Ali) Will an array of this imrpove IO efficiency?

struct Chrom
{
    string chrName;
    long long chrStart;
};

vector <Chrom> chrom;

struct QItem {
    int count[2]; 		// Number of reads covering the base in [0] methylated, and [1] unmethylated form.
    long long pos; 		// The genomic position of one Cytosine in either + or - strand (starting from 1)
    int chr; 			// Chromosome number of each cytosine in the circular queue
} queue[BUFFER_SIZE];

int ChromIndex(char * chr) {
    for(unsigned int i = 0; i < chrom.size(); i++) {
        if(chrom[i].chrName == chr)
            return i;
    }
    return -1;

}

int checkGAorCT() {
    int GA = 0, CT = 0;
    long long pos = chrom[ChromIndex(line.chr)].chrStart + line.pos-1;
    for (unsigned int i = 0; i < line.seq_string.size(); i++)
        if (toupper(line.seq_string[i])=='A' && toupper(genome[pos+i])=='G') GA++;
        else if (toupper(line.seq_string[i])=='T' && toupper(genome[pos+i])=='C') CT++;
    return (GA > CT) ? 1 : ((GA < CT) ? 0 : -1);
}

// Prints the information of one base stored in the queue[index].

void PrintOutput(int index) {
    long long pos = (signed) queue[index].pos;                       // Location in the whole genome
    if (pos < 0) return;
    int chrNum = queue[index].chr;
    string chrName = chrom[chrNum].chrName;
    long long chrPos = pos - chrom[chrNum].chrStart;	// Location in the chromosome
    char base = toupper(genome[pos - 1]);	// Main base
    char next;
    if (base == 'C') next = ((signed) gs > pos) ? toupper(genome[pos]) : '$';
    else if (base == 'G') next = (pos > 1) ? toupper(genome[pos - 2]) : '$';
    else fprintf(stderr, "Bug: The base at genomic locationd %lld is not either C or G\n", pos);
    char strand = (base == 'C') ? '+' : '-';
    bool cpg = (base == 'C' && next == 'G') || (base == 'G' && next == 'C');
    char context[] = "CG";
    if (! cpg) context[1] = 'H';
    int total = queue[index].count[0] + queue[index].count[1];
    int meth = queue[index].count[0];
    if (allCytosines || cpg || meth > 0)
        fprintf(stdout, "%s\t%lld\t%c\t%d\t%d\t%.4f\t%s\n", chrName.c_str(), chrPos, strand,total, meth, (double) meth / total, context);
}

void ProcessMethylation() {
    char methyl='C', unmethyl='T';
    int chrNum = ChromIndex(line.chr);
    long long refPos = chrom[chrNum].chrStart + line.pos-1;

    if (line.strand == '-') {
        methyl='G';
        unmethyl='A';
    }
    for (unsigned int i=0; i < line.seq_string.size() ; i++)
        if(toupper(genome[refPos+i]) == methyl) { // It's a cytosine either in + or - strands
            long long pos = refPos+i;
            int index = (pos)%BUFFER_SIZE;
            if(queue[index].pos == -1 || queue[index].chr != chrNum || queue[index].pos != pos + 1) { // A different genomic location, the process of which is already finished
                if (queue[index].pos > -1 && queue[index].count[0] + queue[index].count[1] > 0)
                    PrintOutput(index);
                queue[index].count[0] = 0;
                queue[index].count[1] = 0;
                queue[index].pos = pos +1;
                queue[index].chr = chrNum;
            }
            if (toupper(line.seq_string[i]) == methyl)
                queue[index].count[0]++;
            else if(toupper(line.seq_string[i]) == unmethyl)
                queue[index].count[1]++;
        }
}

// Handles mismatches and inserts in cigar by converting reads
void convertRead() {
    int j = 0;
    for(unsigned int i=0; i< line.seq_string.size(); i++) {
        if (line.cigar2[i] == 'i')
            line.seq_string.erase(j,1);
        else if (line.cigar2[i] == 'd') {
            line.seq_string.insert(j,"M");
            j++;
        }
        else
            j++;
    }
}

// reverse complement reads
void reverseRead() {
    int i,j=0 ;
    char * copy = new char[line.seq_string.size()];
    line.seq_string.copy(copy, line.seq_string.size(),0);
    for (i = line.seq_string.size()-1; i >= 0 ; i--) {
        switch (copy[i]) {
        case 'A':
            line.seq_string[j]='T';
            break;
        case 'T':
            line.seq_string[j]='A';
            break;
        case 'C':
            line.seq_string[j]='G';
            break;
        case 'G':
            line.seq_string[j]='C';
            break;

        default:
            break;
        }
        j++;
    }

}

// Generates an array from Cigar sequence for ease of computation
void convertCigar(char * cigar, char * cigar2) {
    int pos = 0;
    int value = 0, j , index =0;;
    while (1) {
        if (!isdigit(cigar[pos])) {
            if (value > 0)
                for (j = 0; j < value; j++) {
                    cigar2[index] = cigar[pos];
                    index++;
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

// read chunks from samfile and does computation on them
bool paired = false;
uint64_t last_pos;
int ProcessSamFile(FILE * samFile, FILE * ambFile) {
    char buffer[maxSamFileLineLength];
    char rnext[100], pnext[100], seq_string[maxReadLength], quality_string[maxReadLength]; //(Ali) I kept most of the limitations, but better to update them based on defined constants
    int flag;
    uint32_t mapq;
    long long int tlen;
    bool stop = false;
    buffer[0] = 0;
    while (! stop) {
        do {
            if(fgets(buffer, maxSamFileLineLength, samFile) == NULL) stop = true;
            //else if (buffer[0] != '@');
        } while (! stop && buffer[0] == '@'); // End of file, header lines
        if(stop)
            return -1;
        sscanf(buffer, "%s\t%d\t%s\t%llu\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", line.rname, &flag, line.chr, (unsigned long long *) &(line.pos) ,&mapq, line.cigar,rnext,pnext, &tlen,seq_string,quality_string);
        if (flag & 4) continue; // Unmapped read
        if((flag & 0x1) && line.pos==last_pos)//paired read and unmapped
            continue;
        line.seq_string = seq_string;
        convertCigar(line.cigar, line.cigar2);
        if(debug)
            cerr << "cigar2:  "<<line.cigar2<< endl;
        convertRead();
        //if (flag & 16) reverseRead();
        int result = checkGAorCT();
        if(debug)
            cerr << result<<endl;
        if (result == -1) {
            amb++;
            if (ambFile) fprintf(ambFile, "%s\n", buffer);
            continue; // The read can not be decided to match either PCR product or original sequence.
            // (Ali) we should add some input argument that shows whether PCR amplification is used or not. If not, we can ignore PCR product probability in deciding GAorCT
        }
        // Now Original+ reads have flag=0, result=0.  Original-: flag=16, result=1.  PCR+: flag=16, result=0.  PCR-: flag=0, result=1
        line.strand = (result) ? '-' : '+';
        ProcessMethylation();
    }//while
    return 0;
}


// Reads genome file and save chromosomes' start positions
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
    char fLineMain[10000];
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
            cerr << name << endl;
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


int main(int argc, char *argv[]) {
    char *referenceName = NULL, *samName = NULL, *ambName = NULL;
    FILE * samFile = NULL;

    static struct option long_options[] = {
        { "ref", required_argument, NULL, 'r' },
        { "sam", required_argument, NULL, 's' },
        { "all", no_argument, NULL, 'a' },
        { "amb", required_argument, NULL, 'm' },
        { NULL, 0, NULL, 0}
    };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "r:s:am:", long_options, &option_index)) >= 0) {
        switch (c) {
        case 'r':
            referenceName = strdup(optarg);
            break;
        case 's':
            samName = strdup(optarg);
            break;
        case 'a':
            allCytosines = true;
            break;
        case 'm':
            ambName = strdup(optarg);
            break;
        }
    }

    if (! referenceName || ! samName ) {
        fprintf(stderr, "usage: %s <-r reference genome> <-s input SAM file> [-a] [-m file]\n Use -a for printing the information of all cytosines (+/- strands) which are not methylated nor in CpG context.\n", argv[0]);
        fprintf(stderr, "Use -m followed by SAM file name as the output for the ambiguous reads.\n");
        return 1;
    }
    ReadGenome(referenceName);
    fprintf(stderr, "Complete reading genome.\n");
    samFile = fopen(samName, "r");
    if (!samFile) {
        fprintf(stderr, "Could not open SAM file %s\n", samName);
        return 1;
    }

    for(int k = 0 ; k < BUFFER_SIZE ; k++) queue[k].pos = -1;
    FILE * ambFile = NULL;
    if (ambName) ambFile = fopen(ambName, "w");
    // Main process starts here
    ProcessSamFile(samFile, ambFile);
    // Flushing out the remaining information in the queue
    for(int i=0; i<BUFFER_SIZE ; i++)
        if (queue[i].pos > -1 && queue[i].count[0] + queue[i].count[1] > 0)
            PrintOutput(i);
    if (samFile) fclose(samFile);
    if (ambFile) fclose(ambFile);
    fprintf(stderr, "There were %lld ambiguous reads, with equal C->T and G->A conversions.\n", amb);
    fprintf(stderr, "Finished successfully.\n");
}//main


#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

const long maxGenomeSize = 4e9;
const int maxChromosomeNum = 1000;
const int interval_side = 3000;
long chromPos[maxChromosomeNum];
unsigned long gs;
char chromName[maxChromosomeNum][100];
char * genome;
int chromNum, xChromosome = -1, yChromosome = -1;
char * genomeFile, * annotationFile, * outputFileIslandConsidered, * outputFileContextConsidered, * outputFileComplete;

void swap(int &a, int &b) {
    int c = a;
    a = b;
    b = c;
}

struct Gene {
    char name[20], refseq[20];
    long tss, tts, chrom;
    bool strand;
}; 

inline void ToLower(char * s) {
    for (int i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] += 'a' - 'A';
}

inline void ToUpper(char * s) {
    for (int i = 0; i < strlen(s); i++)
        if (s[i] >= 'a' && s[i] <= 'z') s[i] += 'A' - 'a';
                   
}

int ChromIndex(char * chrom) {
    ToLower(chrom);
    if (strstr(chrom, "chr") == chrom) chrom += 3;
    if (strstr(chrom, "ch") == chrom) chrom += 2;
	if (chrom[0] >= '0' && chrom[0] <= '9' && strlen(chrom) < 3)
		return atoi(chrom) - 1;
    if (strlen(chrom) > 1) return -1;
	if (chrom[0] == 'x') return xChromosome;
	else 
	    if (chrom[0] == 'y') return yChromosome;
		else return -1;
}


void ReadGenome(char * genomeFile) {
	cerr << "Allocating memory..." << endl;
	genome = new char[maxGenomeSize];
	gs = 0;
	chromNum = 0;
	cerr << "Reading genome..." << endl;
	char fLine[10000];
	FILE * f = fopen(genomeFile, "r");
	if (! f) {
		cerr << "Error: Genome file not found or could not be opened" << endl;
		exit(1);
	}	
    while (! feof(f)) {
		int n = fscanf(f, "%s\n", fLine);
		if (n == EOF) break;
		n = strlen(fLine);
		if (fLine[0] == '>') {
			chromPos[chromNum++] = gs;
			if (chromNum > 1) cerr << " Length: " << chromPos[chromNum - 1] - chromPos[chromNum - 2] <<  endl;
			cerr << fLine;
			strcpy(chromName[chromNum - 1], fLine);
            ToLower(fLine);
            if (xChromosome == 0 && strstr(fLine, "x"))
                xChromosome = chromNum - 1;
            if (yChromosome == 0 && strstr(fLine, "y"))
                yChromosome = chromNum - 1;
		} else {
//            ToUpper(fLine);
			memcpy(genome+gs, fLine, n);
			gs += n;
		}
	}
    chromPos[chromNum] = gs;
    cerr << " Length: " << chromPos[chromNum] - chromPos[chromNum - 1] <<  endl;
	fclose(f);
}

// Returns true if the sequence is in valid position and containing no NAs
bool GetSequence(long chr, bool strand, long wStart, long wEnd, char * seq) {
    //pos--; // Converting 1 base position to 0 base 
    bool result = true;
    memset(seq, 'N', wEnd - wStart + 1);
    seq[wEnd - wStart + 1] = 0;    
    long ps, pt, ws;
    ps = wStart - 1;
    pt = wEnd - 1;
    ws = 0;
    if (ps < 0) {
        ws -= ps;
        ps -= ps; // = 0
        result = false;
    }
    if (pt + chromPos[chr] >= chromPos[chr + 1]) { 
        pt = chromPos[chr + 1] - chromPos[chr] - 1;
        result = false;
    }
                          
    memcpy(seq + ws, genome + ps + chromPos[chr], pt - ps + 1);
    if (! strand) {
        char tmp[20000];
        memcpy(tmp, seq, wEnd - wStart + 2);
        long j = 0;
        for (long i = wEnd - wStart; i >= 0; i--)
            switch(tmp[i]) {
                case 'A': case 'a': seq[j++] = 'T'; break;
                case 'C': case 'c': seq[j++] = 'G'; break;
                case 'G': case 'g': seq[j++] = 'C'; break;
                case 'T': case 't': seq[j++] = 'A'; break;
                default: seq[j++] = 'N'; break;
            }
    }
    for (int i = 0; i < strlen(seq); i++) 
        if (seq[i] == 'N' || seq[i] == 'n') { result = false; break; }
    return result;
}

// Converts all CpGs of a sequence in chromosome "chr" at position "wStart:wEnd" to TpG 
// Returns true if the sequence is in valid position and containing no NAs
bool ConvertSequence(long chr, long wStart, long wEnd) {
    //pos--; // Converting 1 base position to 0 base
    bool result = true;
    if (wStart < 1 || wStart + chromPos[chr] - 1 >= chromPos[chr + 1]) return false;

	//cout << chr << "\t" << wStart << "\t" << wEnd << endl; 
	for (unsigned long i = wStart + chromPos[chr] - 1; i < wEnd + chromPos[chr]; i++) {
		//cout << genome[i];
		if ((genome[i] == 'c' || genome[i] == 'C') && (genome[i+1] == 'g' || genome[i+1] == 'G')) // A CpG
			genome[i] = 'T'; // Simulating bisulfite treatment for non-CpG cytosines
    }
	
	//cout << endl;
	return true;
}

void ConvertWholeGenome() {
	for (unsigned long i = 0; i < gs; i++)
        if ((genome[i] == 'c' || genome[i] == 'C') && genome[i+1] != 'g' && genome[i+1] != 'G') // A CpG
            genome[i] = 'T'; // Simulating bisulfite treatment, converting unmethylated Cytosine to Thymine		
}

void ConvertAllC() {
	for (unsigned long i = 0; i < gs; i++)
		if (genome[i] == 'c' || genome[i] == 'C')
			genome[i] = 'T'; //Convert all the Cs to T
}


void ProcessCpGIslands(char * annotationFile) {
    cerr << "Processing CpG island locations from file: " <<  annotationFile << endl;
    ifstream f(annotationFile);
    if (! f.is_open()) {
        cerr << "Error: CpG island locations file not found or could not be opened" << endl;
        exit(1);
    }
    char fLine[10000], chrom[10], strand[10];
    long wStart, wEnd, chr;
    int bin;
   	f.getline(fLine, sizeof(fLine)); // First row
    while (!f.eof()) {
		fLine[0] = 0;
        f.getline(fLine, sizeof(fLine));
		if (! fLine[0]) continue;
		// cerr << fLine << endl;
		wStart = 0;
        sscanf(fLine, "%s %ld %ld", chrom, &wStart, &wEnd);
		if (! wStart) continue;
		// cerr << chrom << '\t' << wStart << '\t' << wEnd << endl;
        chr = ChromIndex(chrom);
		ConvertSequence(chr, wStart, wEnd);
    }
}

void WriteGenome(char * outputFile)
{
	FILE * f = fopen(outputFile, "w");
	char s[100], tmp;
	for (int i = 0; i < chromNum; i++) {
		fprintf(f, "%s\n", chromName[i]);
		unsigned long j;
		for (j = chromPos[i]; j + 50 < chromPos[i + 1]; j+= 50) {
			tmp = genome[j + 50];
			genome[j + 50] = 0;
			fprintf(f, "%s\n", genome + j);
			genome[j + 50] = tmp;
		}
		tmp = genome[chromPos[i + 1]];
		genome[chromPos[i+1]] = 0;
		fprintf(f, "%s\n", genome + j);
		genome[chromPos[i+1]] = tmp;
	}	
	fclose(f);
}

int main(int argc, char * argv[]) {
// Parsing Arugments
	if (argc < 6) { 
        	cerr << "Usage is -g <reference genome> -a <position of CpG islands file> -o <output bisulfite treated genome>" << endl;
        	exit(1);
    }  
    for (int i = 1; i < argc; i++) 
		if (i + 1 != argc) {
            if (strcmp(argv[i], "-g") == 0) 
                genomeFile = argv[++i];
            else 
                if (strcmp(argv[i], "-a") == 0) 
                    annotationFile = argv[++i];
				else 
					if (strcmp(argv[i], "-o") == 0){
						//outputFileIslandConsidered = argv[++i];
						//strcat(outputFileIslandConsidered, "Is")
						outputFileIslandConsidered = "BisulfiteGenomeIslandConsidered";
						outputFileComplete = "BisulfiteGenomeComplete";
						outputFileContextConsidered = "BisulfiteGenomeContextConsidered"
					}
	                else {
    	                cerr << "Not enough or invalid arguments."<< endl;
        	            exit(1);
            	    }
        }
// Reading input
    ReadGenome(genomeFile);
	ConvertWholeGenome();
	WriteGenome(outputFileContextConsidered);
   	ProcessCpGIslands(annotationFile);
	WriteGenome(outputFileIslandConsidered);
	ConvertAllC();
	WriteGenome(outputFileComplete);
	delete [] genome;
}

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <vector>
using namespace std;

vector <long long> chromPos, chromLen;
unsigned long gs;
map <string, int> chromIndex;
map <int, string> chromName, fullName;
// char chromName[maxChromosomeNum][100];
char * genome;
int chromNum, xChromosome = -1, yChromosome = -1;
char * genomeFile, * annotationFile, * CToutputFileIslandConsidered, * CToutputFileContextConsidered, * CToutputFileComplete;
char *originalGenome, * GAoutputFileIslandConsidered, * GAoutputFileContextConsidered, * GAoutputFileComplete;
void swap(int &a, int &b) {
    int c = a;
    a = b;
    b = c;
}

struct window {
    int chr;
    long long wStart, wEnd;
    window(int c, long long s, long long e) {
        chr = c;
        wStart = s;
        wEnd = e;
    }
};

vector <window> islands;

void clean() {
    chromNum = 0;
}

struct Gene {
    char name[20], refseq[20];
    long tss, tts, chrom;
    bool strand;
};

inline void ToLower(char * s) {
    for (unsigned int i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z') s[i] += 'a' - 'A';
}

inline void ToUpper(char * s) {
    for (unsigned int i = 0; i < strlen(s); i++)
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
    else if (chrom[0] == 'y') return yChromosome;
    else return -1;
}

void ReadGenome(string genomeFile) {
    cerr << "Allocating memory..." << endl;
    ifstream ifile(genomeFile.c_str());
    ifile.seekg(0, std::ios_base::end); //seek to end
    //now get current position as length of file
    long long size = ifile.tellg();
    ifile.close();
    genome = new char[size];
    gs = 0;
    chromNum = 0;
    cerr << "Reading genome..." << endl;
    char fLineMain[10000];
    FILE * f = fopen(genomeFile.c_str(), "r");
    if (! f) {
        cerr << "Error: Genome file not found or could not be opened" << endl;
        exit(1);
    }
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
            chromPos.push_back(gs);
            if (chromNum > 0) {
                chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
                cerr << " Length: " << chromLen[chromNum - 1] <<  endl;
            }
            string name = fLine;
			if (name.find(" ") != string::npos) name = name.substr(1, name.find(" ")-1);
			else name = name.substr(1, name.size() - 1);
			cerr << name;
            chromIndex[name] = chromNum;
            chromName[chromNum] = name;
			fullName[chromNum] = fLine;
            chromNum++;
        } else {
            memcpy(genome+gs, fLine, n);
            gs += n;
        }
    }
    chromPos.push_back(gs);
    chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
//  for (int i = 0; i < chromNum; i++)
//      cerr << "ChromPos: " << chromPos[i] << " ChromLen: " << chromLen[i] << endl;
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
            case 'A':
            case 'a':
                seq[j++] = 'T';
                break;
            case 'C':
            case 'c':
                seq[j++] = 'G';
                break;
            case 'G':
            case 'g':
                seq[j++] = 'C';
                break;
            case 'T':
            case 't':
                seq[j++] = 'A';
                break;
            default:
                seq[j++] = 'N';
                break;
            }
    }
    for (unsigned int i = 0; i < strlen(seq); i++)
        if (seq[i] == 'N' || seq[i] == 'n') {
            result = false;
            break;
        }
    return result;
}

// Converts all CpGs of a sequence in chromosome "chr" at position "wStart:wEnd" to TpG
// Returns true if the sequence is in valid position and containing no NAs
bool ConvertSequence(long chr, long wStart, long wEnd, int CT) {
    //pos--; // Converting 1 base position to 0 base
    if (wStart < 1 || wStart + chromPos[chr] - 1 >= chromPos[chr + 1]) {
		cerr << "Warning: the location of CpG island [" << wStart << ", " << wEnd << "] does not match the length of the chromosome " << chr << endl;
		return false;
	}
    //cout << chr << "\t" << wStart << "\t" << wEnd << endl;
    for (unsigned long i = wStart + chromPos[chr] - 1; i < (unsigned) (wEnd + chromPos[chr]); i++) {
        //cout << genome[i];
        if ((genome[i] == 'c' || genome[i] == 'C') && (genome[i+1] == 'g' || genome[i+1] == 'G')) {// A CpG
            if(CT)
                genome[i] = 'T'; // Simulating bisulfite treatment for non-CpG cytosines
            else
                genome[i+1] = 'A';
        }
    }

    //cout << endl;
    return true;
}


// If CT==1: For all C[!G] converts C->T
// If CT==0: For all [!C]G converts G->A
void ConvertWholeGenome(int CT) {
    if(CT) {
        for (unsigned long i = 0; i < gs; i++) {
            if ((genome[i] == 'c' || genome[i] == 'C') && genome[i+1] != 'g' && genome[i+1] != 'G') { // A CpG
                genome[i] = 'T'; // Simulating bisulfite treatment, converting unmethylated Cytosine to Thymine
            }
        }
    } else {
        for (unsigned long i = 0; i < gs; i++) {
            if ((genome[i] != 'c' && genome[i] != 'C') && (genome[i+1] == 'g' || genome[i+1] == 'G')) { // A CpG
                genome[i+1] = 'A'; // Simulating bisulfite treatment, converting unmethylated Guanin to Athenin
            }
        }
    }
}

void ConvertAll(int CT) {
    for (unsigned long i = 0; i < gs; i++) {
        if (CT && (genome[i] == 'c' || genome[i] == 'C'))
            genome[i] = 'T'; //Convert all the Cs to T
        if (!CT && (genome[i] == 'g' || genome[i] == 'G'))
            genome[i] = 'A';
    }
}


void ProcessCpGIslands(char * annotationFile, int CT) {
	bool failed = false;
    cerr << "Processing CpG island locations from file: " <<  annotationFile << endl;
    ifstream f(annotationFile);
    if (! f.is_open()) {
        cerr << "Error: CpG island locations file not found or could not be opened" << endl;
        exit(1);
    }
    char fLine[10000], chrom[10];
    long long wStart, wEnd, chr;
	string last;
    f.getline(fLine, sizeof(fLine)); // First row
    while (!f.eof()) {
        fLine[0] = 0;
        f.getline(fLine, sizeof(fLine));
        if (! fLine[0]) continue;
        // cerr << fLine << endl;
        wStart = 0;
        sscanf(fLine, "%s %lld %lld", chrom, &wStart, &wEnd);
        if (! wStart) continue;
        // cerr << chrom << '\t' << wStart << '\t' << wEnd << endl;
		if (chromIndex.find(chrom) == chromIndex.end()) {
			if (last != chrom) cerr << "Error: chromosome name " << chrom << " does not exist in the reference genome.\nCheck the reference genome FASTA file and make sure the chromosome names match the CpG-islands file.\n";
			last = chrom;
			failed = true;
		}
        chr = chromIndex[chrom];
        islands.push_back(window(chr, wStart, wEnd));
        ConvertSequence(chr, wStart, wEnd, CT);
    }
	if (failed) exit(1);
}

void WriteGenome(char * outputFile)
{
    FILE * f = fopen(outputFile, "w");
	if (! f) {
		fprintf(stderr, "Error opening the output file: %s\n", outputFile);
		exit(1);
	}
    char tmp;
    for (int i = 0; i < chromNum; i++) {
        fprintf(f, "%s\n", fullName[i].c_str());
        unsigned long j;
        for (j = chromPos[i]; j + 50 < (unsigned) chromPos[i + 1]; j+= 50) {
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
        cerr << "Usage is -g <reference genome> -a <position of CpG islands file> -o <output folder>" << endl;
        exit(1);
    }
    for (int i = 1; i < argc; i++)
        if (i + 1 != argc) {
            if (strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
            else if (strcmp(argv[i], "-a") == 0)
                annotationFile = argv[++i];
            else if (strcmp(argv[i], "-o") == 0) {
                char * outputFolder = (char*) malloc(strlen(argv[i+1])+50);
                strcpy(outputFolder, argv[++i]);
                strcat(outputFolder, "/");
                originalGenome = (char*) malloc(strlen(outputFolder)+50);
                strcpy(originalGenome, outputFolder);
                strcat(originalGenome, "originalGenome.fa");
                CToutputFileIslandConsidered = (char*)malloc(strlen(outputFolder)+50);
                strcpy(CToutputFileIslandConsidered, outputFolder);
                strcat(CToutputFileIslandConsidered, "BisulfiteGenomeIslandConsideredCT.fa");
                // CToutputFileContextConsidered = (char*)malloc(strlen(outputFolder)+50);
                // strcpy(CToutputFileContextConsidered, outputFolder);
                // strcat(CToutputFileContextConsidered, "BisulfiteGenomeContextConsideredCT.fa");
                CToutputFileComplete = (char*)malloc(strlen(outputFolder)+50);
                strcpy(CToutputFileComplete, outputFolder);
                strcat(CToutputFileComplete, "BisulfiteGenomeCompleteCT.fa");
                GAoutputFileIslandConsidered = (char*)malloc(strlen(outputFolder)+50);
                strcpy(GAoutputFileIslandConsidered, outputFolder);
                strcat(GAoutputFileIslandConsidered, "BisulfiteGenomeIslandConsideredGA.fa");
                // GAoutputFileContextConsidered = (char*)malloc(strlen(outputFolder)+50);
                // strcpy(GAoutputFileContextConsidered, outputFolder);
                // strcat(GAoutputFileContextConsidered, "BisulfiteGenomeContextConsideredGA.fa");
                GAoutputFileComplete = (char*)malloc(strlen(outputFolder)+50);
                strcpy(GAoutputFileComplete, outputFolder);
                strcat(GAoutputFileComplete, "BisulfiteGenomeCompleteGA.fa");
            }
            else {
                cerr << "Not enough or invalid arguments."<< endl;
                exit(1);
            }
        }
// Reading input
    ReadGenome(genomeFile);
    WriteGenome(originalGenome);
    ConvertWholeGenome(1);
    // WriteGenome(CToutputFileContextConsidered);
    ProcessCpGIslands(annotationFile, 1);
    WriteGenome(CToutputFileIslandConsidered);
    ConvertAll(1);
    WriteGenome(CToutputFileComplete);
    delete [] genome;
    clean();
//Convert G to A
    ReadGenome(genomeFile);
    ConvertWholeGenome(0);
    // WriteGenome(GAoutputFileContextConsidered);
    ProcessCpGIslands(annotationFile, 0);
    WriteGenome(GAoutputFileIslandConsidered);
    ConvertAll(0);
    WriteGenome(GAoutputFileComplete);
    delete [] genome;
}

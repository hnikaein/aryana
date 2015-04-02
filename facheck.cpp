#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

const int maxChromosomeNum = 1000;
const int maxReadLen = 10000;
const int maxSamLineLength = 3 * maxReadLen;
long long chromPos[maxChromosomeNum], chromLen[maxChromosomeNum];
long long gs;
map <string, int> chromIndex;
map <int, string> chromName;
char * genome;
int chromNum = 0;

void ReadGenome(string genomeFile, char * &genome, long long & gs) {
    // cerr << "Allocating memory..." << endl;
    ifstream ifile(genomeFile.c_str());
    ifile.seekg(0, std::ios_base::end);	//seek to end
    //now get current position as length of file
    long long size = ifile.tellg();
    ifile.close();
    genome = new char[size];
    gs = 0;
    chromNum = 0;
    //cerr << "Reading genome..." << endl;
    char fLine[10000];
    FILE * f = fopen(genomeFile.c_str(), "r");
    if (! f) {
        cerr << "Error: Genome file not found or could not be opened" << endl;
        exit(1);
    }
    while (! feof(f)) {
        fLine[0] = 0;
        int n = fscanf(f, "%s\n", fLine);
        if (n == EOF) break;
        n = strlen(fLine);
        if (n == 0) break;
        if (fLine[0] == '>') {
            chromPos[chromNum] = gs;
            if (chromNum > 0) {
                chromLen[chromNum - 1] = chromPos[chromNum] - chromPos[chromNum - 1];
                //		cerr << " Length: " << chromLen[chromNum - 1] <<  endl;
            }

            string name = fLine;
            if (name.find("|") != string::npos) name = name.substr(1, name.find("|")-1);
            else name = name.substr(1, name.size() - 1);
            //	cerr << name;
            chromIndex[name] = chromNum;
            chromName[chromNum] = name;
            chromNum++;
        } else {
            memcpy(genome+gs, fLine, n);
            gs += n;
        }
    }
    chromPos[chromNum] = gs;
    chromLen[chromNum - 1] = chromPos[chromNum] - chromPos[chromNum - 1];
//	for (int i = 0; i < chromNum; i++)
//		cerr << "ChromPos: " << chromPos[i] << " ChromLen: " << chromLen[i] << endl;
//    cerr << " Length: " << chromPos[chromNum] - chromPos[chromNum - 1] <<  endl;
    fclose(f);
}


void revcomp(char * a, long long l = 0) {
    if (! l) l = strlen(a);
    char * b = new char[l];
    memcpy(b, a, l);
    for (long long i = 0; i < l; i++)
        switch (b[i]) {
        case 'a':
        case 'A':
            a[l - i - 1] = 'T';
            break;
        case 'c':
        case 'C':
            a[l - i - 1] = 'G';
            break;
        case 'g':
        case 'G':
            a[l - i - 1] = 'C';
            break;
        case 't':
        case 'T':
            a[l - i - 1] = 'A';
            break;
        default:
            a[l-i-1] = 'N';
        };
    delete [] b;
}


void GetSequence(long long start, long long end, char * seq) {
    if (end < gs) memcpy(seq, genome + start, end - start + 1);
    else memcpy(seq, genome + 2*gs - end - 1, end - start + 1);
    seq[end - start + 1] = 0;
    if (end >= gs) revcomp(seq, end - start + 1);
}

int main(int argc, char * argv[]) {
    long long i, p, pos[2], l, err=0;
    char * seq[2];
    if (argc < 3) {
        cerr << "Usage: facheck <reference genome> <aryana genome>" << endl;
        exit(1);
    }

    ReadGenome(argv[1], genome, gs);
    char * genome2;
    long long gs2;
    ReadGenome(argv[2], genome2, gs2);
//	if (gs2 != 2 * gs)
    cerr << "Genome lengths, gs=" << gs << ", gs2=" << gs2 << endl;
    for (long long i = 0; i < gs; i++)
        if (genome[i] != genome2[i]) cerr << "Pos " << i << " " << genome[i] << "->" << genome2[i] << endl;
    revcomp(genome);
    for (long long i = 0; i < gs; i++)
        if (genome[i] != genome2[gs+i]) cerr << "Pos " << gs+i << " " << genome[i] << "->" << genome2[i] << endl;
    return 0;
}

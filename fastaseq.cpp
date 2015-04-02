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
string genomeFile;

void ReadGenome(string genomeFile) {
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

// Start is 1-based.
bool GetSequence(string chr, long long start, long long length, char *seq) {
    start--;
    if (chromIndex.find(chr) == chromIndex.end() || chromLen[chromIndex[chr]] < start + length) return false;
    memcpy(seq, genome + chromPos[chromIndex[chr]] + start, length);
    seq[length] = 0;
    return true;
}

void PrintSequence(string chr, long long start, long long end, bool revComp) {
    char * seq = new char[end - start + 2];
    if (chr != "") {
        if (! GetSequence(chr, start, end, seq)) {
            cerr << "Error in the given chromosomal position." << endl;
            return;
        }
    }
    else {
        //cerr << start << "  " << end << " " << gs << endl;
        if (start < 0 || end >= 2 * gs)
        {
            cerr << "Error in the given chromosomal position." << endl;
            return;
        }
        if (end < gs) memcpy(seq, genome + start, end - start + 1);
        else memcpy(seq, genome + 2*gs - end - 1, end - start + 1);
        seq[end - start + 1] = 0;
    }
    if (revComp || end >= gs) revcomp(seq, end - start + 1);
    cout << seq << endl;
}

int main(int argc, char * argv[]) {
// Parsing Arugments
    bool revComp = false;
    string chr;
    long long start = 0, end = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-r") == 0)
            revComp = true;

        if (i < argc - 1) { // The paired arguments
            if (strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
            else if (strcmp(argv[i], "-sl") == 0 && i + 2 < argc) {
                start = atoi(argv[++i]);
                end = start + atoi(argv[++i]) - 1;
            } else if (strcmp(argv[i], "-se") == 0 &&  i + 2 < argc) {
                start = atoi(argv[++i]);
                end = atoi(argv[++i]);
            } else if (strcmp(argv[i], "-c") == 0 && i + 3 < argc) {
                chr = argv[++i];
                start = atoi(argv[++i]);
                end = atoi(argv[++i]);
            }
            else {
                cerr << "Unrecognized argument or invalid usage: "<< argv[i] << endl;
                exit(1);
            }
            continue;
        }
        cerr << "This argument should be followed by an input: "<< argv[i] << endl;
        exit(1);
    }
    if (genomeFile == "" || ! start) {
        cerr << "Arguments: -g <reference genome, mandatory>" << endl <<
             "-sl start length (absolute 0-based start position)" << endl <<
             "-se start end (absolute 0-based start and end positions)" << endl <<
             "-c chr start end (1-based start and end)\n-r (reverse-complement)" << endl;
        exit(1);
    }

// Reading input
    ReadGenome(genomeFile);
    PrintSequence(chr, start, end, revComp);
    delete [] genome;
}

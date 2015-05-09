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

const int maxReadLen = 10000;
const int maxSamLineLength = 3 * maxReadLen;
vector <long long> chromPos, chromLen;
long long gs;
map <string, int> chromIndex;
map <int, string> chromName;
char * genome;
int chromNum = 0;
string genomeFile;
enum inputType {startLength, startEnd, chrLoc} inpType;

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
                //		cerr << " Length: " << chromLen[chromNum - 1] <<  endl;
            }

            string name = fLine;
            if (name.find(" ") != string::npos) name = name.substr(1, name.find(" ")-1);
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
    chromPos.push_back(gs);
    chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
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
            delete[] seq;
            return;
        }
    }
    else {
        //cerr << start << "  " << end << " " << gs << endl;
        if (start < 0 || end >= 2 * gs)
        {
            cerr << "Error in the given chromosomal position." << endl;
            delete [] seq;
            return;
        }
        if (end < gs) memcpy(seq, genome + start, end - start + 1);
        else memcpy(seq, genome + 2*gs - end - 1, end - start + 1);
        seq[end - start + 1] = 0;
    }
    if (revComp || end >= gs) revcomp(seq, end - start + 1);
    cout << seq << endl;
    delete[] seq;
}

int main(int argc, char * argv[]) {
// Parsing Arugments
    bool revComp = false;
    string chr;
    long long start = 0, end = 0, length;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-r") == 0)
            revComp = true;
        else if (strcmp(argv[i], "-sl") == 0)
            inpType = startLength;
        else if (strcmp(argv[i], "-se") == 0)
            inpType = startEnd;
        else if (strcmp(argv[i], "-c") == 0)
            inpType = chrLoc;
        else if ((i < argc - 1) && strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
        else {
            cerr << "Unrecognized argument or invalid usage: "<< argv[i] << endl;
            exit(1);
        }
    }
    if (genomeFile == "") {
        cerr << "Usage: fastaseq -g <reference genome, mandatory> {-sl, -se, -c} {-r for reverse-complementing output}" << endl <<
             "Genomic positions are then read from the standard input, according to the given argument as below:" << endl <<
             "-sl (start length):         two columns of integers indicating absolute 0-based start position and length of the region" << endl <<
             "-se (start end):            two columns of integers indicating absolute 0-based start and end of the region" << endl <<
             "-c  (chromosomal location): a column of chromosome names followed by two columns of integers indicating 1-based start and end positions inside chromosome" << endl <<
             "If absolute positions are beyong genome size, they are treated as the locations in reverse-complemented genome concatenated to the original genome" << endl;
        exit(1);
    }

// Reading input
    ReadGenome(genomeFile);
    while (cin.good()) {
        start = -1;
        switch (inpType) {
        case startLength:
            cin >> start >> length;
            end = start + length - 1;
            break;
        case startEnd:
            cin >> start >> end;
            break;
        case chrLoc:
            cin >> chr >> start >> end;
            break;
        };
        if (start < 0) break;
        PrintSequence(chr, start, end, revComp);
    }
    delete [] genome;
}

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
    if (argc < 4) {
        cerr << "Usage: checkbwt <reference genome fasta file> <BWT index positions> <length of strings to check> <-a for just listing all sequences>" << endl;
        exit(1);
    }

    ReadGenome(argv[1]);
    l = atoi(argv[3]);
    for (i = 0; i < 2; i++) seq[i] = new char[l + 1];
    seq[1][0] = 0;
    FILE * f = fopen(argv[2], "r");
    bool listall = (argc == 5) && (strcmp(argv[4], "-a") == 0);
    for (i = 0; i < 2 * gs; i++)
    {
        fscanf(f, "%lld", &p);
        pos[i&1] = p;
        if (p < gs && p + l - 1 >= gs) {
            GetSequence(p, gs - 1, seq[i&1]);
            GetSequence(gs, l + p - 1, seq[i&1]+gs-p);
        } else if (p + l > 2 * gs) {
            //cerr << "[" << p << "," << 2*gs-1 << "], [0, " << l-2*gs+p-1<<"]\n";
            GetSequence(p, 2*gs-1, seq[i&1]);
            //GetSequence(0, l- 2 * gs + p - 1, seq[i&1] + 2 * gs - p);
        } else GetSequence(p, p + l - 1, seq[i&1]);
        if (listall) printf("%s\n", seq[i&1]);
        else if (strcmp(seq[i&1],seq[1-(i&1)]) < 0) {
            printf("[index: %lld, pos: %lld, seq: %s] > [index: %lld, pos: %lld, seq: %s]\n", i-1, pos[1-(i&1)], seq[1-(i&1)], i, p, seq[i&1]);
            err++;
        }
    }
    delete [] genome;
    fclose(f);
    if (! listall) cerr << "Identified " << err << " errors in the given BWT table.\n";
    return 0;
}

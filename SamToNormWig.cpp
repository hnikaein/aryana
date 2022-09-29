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

const int maxChromosomeNum = 1000000;
const int maxReadLen = 10000;
const int maxSamLineLength = 3 * maxReadLen;
bool strandSpecific = false, normalize = false, endPlusPlus = true;
int * coveragePos[maxChromosomeNum], * coverageNeg[maxChromosomeNum];
long long chromSize[maxChromosomeNum];
map <string, int> chromIndex;
map <int, string> chromName;
int chromNum = 0;
string chromSizesFile, samFileName, outputFileName, negativeOutputFileName;
double scaleFactor = 1e8, totalCoverage = 0;

void ReadChromSizes(string chromSizesFile) {
    FILE * f = fopen(chromSizesFile.c_str(), "r");
    if (! f) {
        cerr << "Error: Chromosome sizes file not found or could not be opened" << endl;
        exit(1);
    }
    char name[1000];
    long long len;
    while (! feof(f)) {
        len = 0;
        if (!fscanf(f, "%s %lld", name, &len)) break;
        if (! len) break;
        chromIndex[name] = chromNum;
        chromName[chromNum] = name;
        chromSize[chromNum] = len;
        coveragePos[chromNum] = new int[len];
        bzero(coveragePos[chromNum], len * sizeof(int));
        if (strandSpecific) {
            coverageNeg[chromNum] = new int[len];
            bzero(coverageNeg[chromNum], len * sizeof(int));
        }
        cerr << "Chromosome " << name << ", length: " << len << endl;
        chromNum++;
    }
    fclose(f);
}

// Returns the total length of reference sequence pointed to by the given cigar sequence

void ProcessSamRecord(string qname, int flag, string chr, long long pos, string cigar) {
    pos--;
    unsigned int cigarpos = 0;
    int ch = chromIndex[chr];
    long long p;
    while (cigarpos < cigar.size()) {
        int num = 0;
        while (cigarpos < cigar.size() && cigar[cigarpos] >= '0' && cigar[cigarpos] <= '9') num = num * 10 + cigar[cigarpos++] - '0';
        if (! num || cigarpos >= cigar.size()) {
            cerr << "Invalid cigar sequence " << cigar << endl;
            return;
        }
        switch (tolower(cigar[cigarpos++])) {
        case 'm':
            for (p = pos; p < pos + num; p++) {
                if (p >= chromSize[ch]) {
                    cerr << "Error: read " << qname << " mapping position is higher than the length of chromosome " << chr << endl;
                    return;
                }
                else {
                    if (! strandSpecific || (flag & 16) == 0) coveragePos[ch][p]++;
                    else coverageNeg[ch][p]++;
                }
            }
            pos += num;
            totalCoverage += num;
            break;
        case 'd':
            pos += num;
        }
    }
}

void PrintStrand(FILE * f, int ** coverage) {
    for (int i = 0; i < chromNum; i++) {
        const char * ch = chromName[i].c_str();
        long long p = 0, pos = 0, p2;
        while (p < chromSize[i]) {
            while (pos < chromSize[i] && coverage[i][pos] == coverage[i][p]) pos++;
            p2 = (endPlusPlus) ? pos+1 : pos;
            if (p2 > chromSize[i]) p2 = chromSize[i];
            if (coverage[i][p]) {
                if (normalize)
                    fprintf(f, "%s\t%lld\t%lld\t%f\n", ch, p+1, p2, coverage[i][p] * scaleFactor);
                else if (scaleFactor > 0) fprintf(f, "%s\t%lld\t%lld\t%d\n", ch, p+1, pos + 1, coverage[i][p]);
                else fprintf(f, "%s\t%lld\t%lld\t%d\n", ch, p+1, p2, -coverage[i][p]);
            }
            p = pos;
        }
    }
}

void PrintOutput() {
    scaleFactor = scaleFactor / totalCoverage;
    FILE * f = stdout;
    if (outputFileName != "") f = fopen(outputFileName.c_str(), "w");
    PrintStrand(f, coveragePos);
    if (strandSpecific) {
        scaleFactor *= -1;
        if (negativeOutputFileName == "")
            PrintStrand(f, coverageNeg);
        else {
            fclose(f);
            f = fopen(negativeOutputFileName.c_str(), "w");
            PrintStrand(f, coverageNeg);
        }
    }
    fclose(f);
}

void ProcessSamFile(string samFileName) {
    long long int tlen, pos;
    int flag, mapq;
    char qname[1000], rname[1000], rnext[1000], pnext[1000], seq[maxReadLen], quality_string[maxReadLen], cigar[2*maxReadLen];
    //copy[maxReadLen], qname2[1000], rnext2[1000], pnext2[1000], seq2[maxReadLen], quality_string2[maxReadLen], copy2[maxReadLen], cigar[2*maxReadLen], cigar2[2*maxReadLen], refSeq[2*maxReadLen], refSeq2[2*maxReadLen];
    FILE * f = stdin;
    if (samFileName != "") f = fopen(samFileName.c_str(), "r");
    if (! f) {
        cerr << "Error opening SAM file" << endl;
        exit(-1);
    }
    char buf[maxSamLineLength];
    while (! feof(f)) {
        buf[0] = 0;
        if (! fgets(buf, sizeof(buf), f) || !buf[0]) break;
        if (buf[0] == '@') continue;
//		cerr << "SAM line: " << buf << endl;
        sscanf(buf,"%s\t%d\t%s\t%lld\t%d\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos, &mapq, cigar, rnext, pnext, &tlen, seq, quality_string);
//		cerr << "Position: " << rname << ":" << pos << endl;
        if (strcmp(rname, "*") != 0) ProcessSamRecord(qname, flag, rname, pos, cigar);
    }
}


int main(int argc, char * argv[]) {
// Parsing Arugments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0) {
            strandSpecific = true;
            continue;
        }
        if (strcmp(argv[i], "-n") == 0) {
            normalize = true;
            continue;
        }
        if (i < argc - 1) { // The paired arguments
            if (strcmp(argv[i], "-c") == 0)
                chromSizesFile = argv[++i];
            else if (strcmp(argv[i], "-i") == 0)
                samFileName = argv[++i];
            else if (strcmp(argv[i], "-w") == 0)
                outputFileName = argv[++i];
            else if (strcmp(argv[i], "-W") == 0) {
                negativeOutputFileName = argv[++i];
                strandSpecific = true;
            } else if (strcmp(argv[i], "-e")==0)
                endPlusPlus=false;
            else if (strcmp(argv[i], "-f") == 0) {
                normalize = true;
                scaleFactor = atof(argv[++i]);
            }
            else {
                cerr << "Unrecognized argument: "<< argv[i] << endl;
                exit(1);
            }
            continue;
        }
        cerr << "This argument should be followed by an input: "<< argv[i] << endl;
        exit(1);
    }
    if (chromSizesFile == "") {
        cerr << "Arguments: -c <chromosome sizes, mandatory> -i <alignment sam file> -w <result wig file> -s (strand-specific output, the WIG file will have both positive and negative values)" << endl <<
             "           -n (normalize, default=no) -f <scaling factor, default=1e+8, sets -n on> -W <negative strand wig file, sets -s on, without it both strands are written to the -o specified file>" << endl <<
             "           -e (end of each interval will be less than start of the next interval. default: end of each interval can be equal to start of the next interval)" << endl;
        exit(1);
    }

// Reading input
    ReadChromSizes(chromSizesFile);
    cerr << "Processing SAM file..." << endl;
    ProcessSamFile(samFileName);
    PrintOutput();
    cerr << "Finished." << endl;
    chromIndex.clear();
    chromName.clear();
    if (strandSpecific) for (int i = 0; i < chromNum; i++) delete[] coverageNeg[i];
    for (int i = 0; i < chromNum; i++) delete[] coveragePos[i];
}

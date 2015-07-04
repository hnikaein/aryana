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
string genomeFile, samFileName, outputFileName;
bool bisSeq = false, realPos = false, seqOutput = false;

void ReadGenome(string genomeFile) {
    cerr << "Allocating memory..." << endl;
    ifstream ifile(genomeFile.c_str());
    ifile.seekg(0, std::ios_base::end);	//seek to end
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
            chromNum++;
        } else {
            memcpy(genome+gs, fLine, n);
            gs += n;
        }
    }
    chromPos.push_back(gs);
    chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
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

// Start is 1-based.
bool GetSequence(string chr, long long start, long long length, char *seq) {
    start--;
    if (chromIndex.find(chr) == chromIndex.end() || chromLen[chromIndex[chr]] < start + length) return false;
    memcpy(seq, genome + chromPos[chromIndex[chr]] + start, length);
    seq[length] = 0;
    return true;
}

// Returns the total length of reference sequence pointed to by the given cigar sequence

int CigarLen(string cigar) {
    unsigned int cigarpos = 0, len = 0;
    while (cigarpos < cigar.size()) {
        int num = 0;
        while (cigarpos < cigar.size() && cigar[cigarpos] >= '0' && cigar[cigarpos] <= '9') num = num * 10 + cigar[cigarpos++] - '0';
        if (! num || cigarpos >= cigar.size()) {
            cerr << "Invalid cigar sequence " << cigar << endl;
            return -1;
        }
        switch (tolower(cigar[cigarpos++])) {
        case 'm':
        case 'd':
            len += num;
        }
    }
    return len;
}


// Finds the best alignment between a (reference) and b (query).
// Permits conversions from orig_base inside a to conv_base inside b with 0 penalty. returns the edit distance

int EditDistance(string a, string b, string cigar, char orig_base = '\0', char conv_base = '\0') {
    //cerr << "Edit distance, length(a): " << a.size() << ", length(b): " << b.size() << ", cigar: " << cigar << endl << "a: " << a << endl << "b: " << b << endl;
    unsigned int apos = 0, bpos = 0, cigarpos = 0;
    int penalty = 0;
    while (cigarpos < cigar.size()) {
        int num = 0;
        while (cigarpos < cigar.size() && cigar[cigarpos] >= '0' && cigar[cigarpos] <= '9') num = num * 10 + cigar[cigarpos++] - '0';
        if (! num || cigarpos >= cigar.size()) {
            cerr << "Invalid cigar sequence " << cigar << endl;
            return -1;
        }
        switch (tolower(cigar[cigarpos++])) {
        case 'm':
            if (apos + num > a.size() || bpos + num > b.size()) {
                cerr << "Invalid cigar length compared to the sequences, cigar: " << cigar << ", apos+num: " << apos + num << ", a.size: " << a.size() << ", bpos+num: " << bpos + num << ", b.size: " << b.size() << endl;
                return -1;
            }
            for (int i = 0; i < num; i++)
                if (toupper(a[apos+i]) != toupper(b[bpos+i]) && (toupper(a[apos+i]) != orig_base || toupper(b[bpos+i]) != conv_base)) {
                    //cerr << "Mismatch at " << i << ", a[" << apos + i << "]=" << a[apos+i] << ", b[" << bpos+i << "]=" << b[bpos+i] << endl;
                    penalty++;
                }
            apos += num;
            bpos += num;
            break;
        case 'i':
            if (bpos + num > b.size()) {
                cerr << "Invalid cigar length compared to the sequences " << cigar << endl;
                return -1;
            }
            bpos += num;
            penalty += num;
            break;
        case 'd':
            if (apos + num > a.size()) {
                cerr << "Invalid cigar length compared to the sequences " << cigar << endl;
                return -1;
            }
            apos += num;
            penalty += num;
            break;
        }
    }
    //cerr << "Penlty: " << penalty << endl;
    return penalty;
}

void PrintOutput(FILE * f, string name = "NA", string chr1 = "NA", long long pos1 = 0, int penalty1 = 0, string chr2 = "NA", long long pos2 = 0, int penalty2 = 0, char * samSeq = NULL, char * cigar = NULL, char * seq1 = NULL, char* seq2 = NULL) {
    if (! realPos) {
        fprintf(f, "%s\t%s\t%lld\t%d", name.c_str(), chr1.c_str(), pos1, penalty1);
        if (seqOutput) fprintf(f, "\t%s\t%s\t%s", samSeq, cigar, seq1);
    }
    else {
        fprintf(f, "%s\t%s\t%lld\t%d\t%s\t%lld\t%d", name.c_str(), chr1.c_str(), pos1, penalty1, chr2.c_str(), pos2, penalty2);
        if (seqOutput) fprintf(f, "\t%s\t%s\t%s\t%s", samSeq, cigar, seq1, seq2);
    }
    fprintf(f, "\n");
}

void ProcessSamFile(string samFileName) {
    long long int tlen, pos;
    int flag, mapq, penalty, penalty2, startPos;
    char qname[1000], rname[1000], rnext[1000], pnext[1000], seq[maxReadLen], quality_string[maxReadLen], cigar[2*maxReadLen], refSeq[2*maxReadLen], refSeq2[2*maxReadLen];
    FILE * f = stdin, * of = stdout;
    if (samFileName != "") f = fopen(samFileName.c_str(), "r");
    if (outputFileName != "") of = fopen(outputFileName.c_str(), "w");
    if (! f) {
        cerr << "Error opening SAM file" << endl;
        exit(-1);
    }
    char buf[maxSamLineLength];
    if (realPos) {
        fprintf(of, "ReadName\tAlnChr\tAlnPos\tAlnPen\tRealChr\tRealPos\tRealPen");
        if (seqOutput) fprintf(of, "\tSamSeq\tCIGAR\tAlnSeq\tRealSeq");
    }
    else {
        fprintf(of, "ReadName\tAlnChr\tAlnPos\tAlnPen");
        if (seqOutput) fprintf(of, "\tSamSeq\tCIGAR\tAlnSeq");
    }
    fprintf(of, "\n");
    while (! feof(f)) {
        buf[0] = 0;
        penalty = penalty2 = 1000000;
        startPos = 0;
        string tmp, chr = "NA", start, end, type;
        if (! fgets(buf, sizeof(buf), f) || !buf[0]) break;
        if (buf[0] == '@') continue;
//		cerr << "SAM line: " << buf << endl;
        sscanf(buf,"%s\t%d\t%s\t%lld\t%d\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, rname, &pos, &mapq, cigar, rnext, pnext, &tlen, seq, quality_string);
//		cerr << "Position: " << rname << ":" << pos << endl;
        if (strcmp(rname, "*") == 0) {
            strcpy(rname, "NA");
            strcpy(refSeq, "NA");
            pos = 0;
        } else if (GetSequence(rname, pos, CigarLen(cigar), refSeq)) {
//			cerr << "Flag: " << flag << ", Read: " << seq << ", Ref: " << refSeq << endl;
//			if (flag & 16) {
//				revcomp(refSeq);
//				cerr << "Converted Ref: " << refSeq << endl;
//			}
            if (bisSeq)
                penalty = min(EditDistance(refSeq, seq, cigar, 'C', 'T'), EditDistance(refSeq, seq, cigar, 'G', 'A'));
            else
                penalty = EditDistance(refSeq, seq, cigar);
        }
        if (realPos) {
            // BisSimul simulated reads
            stringstream s(qname);
            getline(s, tmp, '|');
            getline(s, chr, ':');
            getline(s, start, '-');
            getline(s, end, '|');
            getline(s, type, ' ');
            if (GetSequence(chr, atoi(start.c_str()), atoi(end.c_str()) - atoi(start.c_str()) + 1, refSeq2)) {
                startPos = atoi(start.c_str());
                if (type == "-o" || type == "+p") {
                    revcomp(refSeq2);
                }
                if (type[1] == 'o') penalty2 = EditDistance(refSeq2, seq, "100m", 'C', 'T');
                else penalty2 = EditDistance(refSeq2, seq, "100m", 'G', 'A');
            }
        }
        PrintOutput(of, qname, rname, pos, penalty, chr, startPos, penalty2, seq, cigar, refSeq, refSeq2);
    }
}


int main(int argc, char * argv[]) {
// Parsing Arugments
    for (int i = 1; i < argc; i++) {
        if (i < argc - 1) { // The paired arguments
            if (strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
            else if (strcmp(argv[i], "-i") == 0)
                samFileName = argv[++i];
            else if (strcmp(argv[i], "-o") == 0)
                outputFileName = argv[++i];
            else {
                cerr << "Unrecognized argument: "<< argv[i] << endl;
                exit(1);
            }
            continue;
        }
        cerr << "This argument should be followed by an input: "<< argv[i] << endl;
        exit(1);
    }
    if (genomeFile == "") {
        cerr << "Arguments: -g <reference genome, mandatory> -i <alignment sam file> -o <result file> -S (print SAM, CIGAR and alignment sequences in output)" << endl \
             << "-b (bisulfite-seq input, default: DNA-seq)  -s (use only if reads are generated via BisSimul or contain the true positions in a similar way)" << endl;
        exit(1);
    }

// Reading input
    cerr << "Reading reference genome..." << endl;
    ReadGenome(genomeFile);
    cerr << "Processing SAM file..." << endl;
    ProcessSamFile(samFileName);
    cerr << "Finished." << endl;
    delete [] genome;
}

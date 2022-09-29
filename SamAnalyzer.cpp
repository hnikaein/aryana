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

vector <long long> chromPos, chromLen;
long long gs;
map <string, int> chromIndex;
map <int, string> chromName;
char * genome;
int chromNum = 0;
string genomeFile, samFileName, outputFileName;
bool bisSeq = false, realPos = false, seqOutput = false, onlyPenalties = false, keepSecondary = false;

void ReadGenome(string genomeFile) {
    cerr << "Allocating memory..." << endl;
    ifstream ifile(genomeFile.c_str());
    ifile.seekg(0, std::ios_base::end);	//seek to end
    //now get current position as length of file
    long long size = ifile.tellg();
    ifile.close();
	if (size <= 0) {
		cerr << "Error reading the genome file: " << genomeFile << endl;
		exit(1);
	}
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
                //cerr << " Length: " << chromLen[chromNum - 1] <<  endl;
            }

            string name = fLine;
            if (name.find(" ") != string::npos) name = name.substr(1, name.find(" ")-1);
            else name = name.substr(1, name.size() - 1);
            //cerr << name;
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
    //cerr << " Length: " << chromPos[chromNum] - chromPos[chromNum - 1] <<  endl;
    fclose(f);
}

void revcomp(string & a, long long l = 0) {
    if (! l) l = a.size();
    string b = a;
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
}

// Start is 1-based.
bool GetSequence(string chr, long long start, long long length, string & seq) {
    start--;
    if (chromIndex.find(chr) == chromIndex.end() || chromLen[chromIndex[chr]] < start + length) return false;
    seq = string(genome + chromPos[chromIndex[chr]] + start, genome + chromPos[chromIndex[chr]] + start + length);
    return true;
}

// Returns the total length of reference sequence pointed to by the given cigar sequence

int CigarLen(string cigar, int & gapOpen, int & gapExt) {
    unsigned int cigarpos = 0, len = 0;
    gapOpen = gapExt = 0;
    while (cigarpos < cigar.size()) {
        int num = 0;
        while (cigarpos < cigar.size() && cigar[cigarpos] >= '0' && cigar[cigarpos] <= '9') num = num * 10 + cigar[cigarpos++] - '0';
        if (! num || cigarpos >= cigar.size()) {
            cerr << "Invalid cigar sequence " << cigar << endl;
            return -1;
        }
        switch (tolower(cigar[cigarpos++])) {
        case 'm':
            len += num;
            break;
        case 'd':
            len += num;
            gapOpen++;
            gapExt += num - 1;
            break;
        case 'i':
        case 's':
        case 'h':
            gapOpen++;
            gapExt += num - 1;
            break;
        }
    }
    return len;
}

string reverse_cigar(string cigar)
{
    int i=0, cigar_size = cigar.size();
    char cigar_tmp[cigar_size + 50];
    int num=0;
    char type=cigar[cigar_size-1];
    int offset=1;
    int lasttmpsize=0;
    for (i=cigar_size-2; i>=0; i--)
    {
        if (cigar[i]>='0' && cigar[i]<='9' )
        {
            num+=offset*(cigar[i]-'0');
            offset*=10;
        }
        else
        {
            lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
            type=cigar[i];
            num=0;
            offset=1;
        }
    }
    lasttmpsize += snprintf(cigar_tmp + lasttmpsize, 10, "%d%c", num,type);
    return cigar_tmp;
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
        case 'h':
        case 's':
            if (bpos + num > b.size()) {
                cerr << "Invalid cigar length compared to the sequences " << cigar << endl;
                return -1;
            }
            bpos += num;
            break;
        case 'd':
            if (apos + num > a.size()) {
                cerr << "Invalid cigar length compared to the sequences " << cigar << endl;
                return -1;
            }
            apos += num;
            break;
        }
    }
    //cerr << "Penlty: " << penalty << endl;
    return penalty;
}

string upper(string s) {
	for (unsigned int i = 0; i < s.size(); i++) s[i] = toupper(s[i]);
	return s;
}

void PrintOutput(FILE * f, string name = "NA", bool aligned = false, int refLength = 0, string chr1 = "NA", long long pos1 = 0, int mismatch1 = 0, int gapOpen1 = 0, int gapExt1 = 0, string chr2 = "NA", long long pos2 = 0, int mismatch2 = 0, int gapOpen2 = 0, int gapExt2 = 0, string samSeq = NULL, string cigar = NULL, string seq1 = NULL, string cigar2 = NULL, string seq2 = NULL) {
    bool matchTogether = realPos && chr1 == chr2 && pos1 == pos2;
    if (onlyPenalties) {
        fprintf(f, "%d\t%d\t%d\t%d\t%d", aligned, refLength, mismatch1, gapOpen1, gapExt1);
        if (realPos) fprintf(f, "\t%d\t%d\t%d\t%d", matchTogether, mismatch2, gapOpen2, gapExt2);
    } else {
        fprintf(f, "%s\t%d\t%s\t%lld\t%d\t%d\t%d\t%d", name.c_str(), aligned, chr1.c_str(), pos1, refLength, mismatch1, gapOpen1, gapExt1);
        if (realPos) fprintf(f, "\t%s\t%lld\t%d\t%d\t%d\t%d", chr2.c_str(), pos2, matchTogether, mismatch2, gapOpen2, gapExt2);
        if (seqOutput) {
            fprintf(f, "\t%s\t%s\t%s", upper(samSeq).c_str(), cigar.c_str(), upper(seq1).c_str());
            if (realPos) fprintf(f, "\t%s\t%s", cigar2.c_str(), upper(seq2).c_str());
        }
    }
    fprintf(f, "\n");
}


void ProcessSamFile(string samFileName) {
    long long int tlen, pos;
    int flag, mapq, mismatch, mismatch2, gapOpen, gapOpen2, gapExt, gapExt2, startPos;
    string qname, rname, rnext, pnext, seq, quality_string, cigar, refSeq, refSeq2, buf;
    if (samFileName != "") 
	if (! freopen(samFileName.c_str(), "r", stdin)) {
        cerr << "Error opening SAM file" << endl;
        exit(-1);
    }
    if (outputFileName != "") 
	if (! freopen(outputFileName.c_str(), "w", stdout)) {
        cerr << "Error opening output file" << endl;
        exit(-1);
    }
    if (onlyPenalties) {
        cout << "Aligned\tRefLength\tAlnMismatch\tAlnGapOpen\tAlnGapExt";
        if (realPos) cout << "\tCorrectPos\tRealMismatch\tRealGapOpen\tRealGapExt";
    } else {
        cout << "ReadName\tAligned\tAlnChr\tAlnPos\tRefLength\tAlnMismatch\tAlnGapOpen\tAlnGapExt";
        if (realPos) cout << "\tRealChr\tRealPos\tCorrectPos\tRealMismatch\tRealGapOpen\tRealGapExt";
        if (seqOutput) {
            cout << "\tSamSeq\tAlnCIGAR\tAlnSeq";
            if (realPos) cout << "\tRealCIGAR\tRealSeq";
        }
    }
    cout << "\n";
    while (true) {
	int cigLen = 0;
        bool aligned = true;
        mismatch = mismatch2 = gapOpen = gapOpen2 = gapExt = gapExt2 = 0;
        startPos = 0;
        string tmp, chr = "NA", start, end, type, cigar2;
	getline(cin, buf);
        if (buf.empty()) break;
        if (buf[0] == '@') continue;
//		cerr << "SAM line: " << buf << endl;
	istringstream bufstream(buf);
        bufstream >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> quality_string;
//		cerr << "Position: " << rname << ":" << pos << endl;
        if (! keepSecondary && (flag & (2048 + 256))) continue;
        if ((flag & 4) || cigar == "*" || rname == "*") {
	    rname = "NA"; 
	    refSeq = "NA";
            pos = 0;
            aligned = false;
        } else {
            cigLen = CigarLen(cigar, gapOpen, gapExt);
            if (cigLen <= 0) continue;
            if (! GetSequence(rname, pos, cigLen, refSeq)) continue;
            if (bisSeq)
                mismatch = min(EditDistance(refSeq, seq, cigar, 'C', 'T'), EditDistance(refSeq, seq, cigar, 'G', 'A'));
            else
                mismatch = EditDistance(refSeq, seq, cigar);
            if (mismatch < 0) continue;
        }
        if (realPos) {
            stringstream s(qname);
            getline(s, tmp, '|');
            getline(s, chr, ':');
            getline(s, start, '-');
            getline(s, end, '|');
            getline(s, type, '|');
            getline(s, cigar2, ' ');
            startPos = atoi(start.c_str());
            if (GetSequence(chr, startPos, CigarLen(cigar2, gapOpen2, gapExt2), refSeq2)) {
                bool samRev = flag & 16, refRev = type == "-o" || type == "+p";
                if (samRev != refRev) revcomp(refSeq2);
                if (samRev && refRev) cigar2 = reverse_cigar(cigar2);
                if (bisSeq)
                    mismatch2 = min(EditDistance(refSeq2, seq, cigar2.c_str(), 'C', 'T'), EditDistance(refSeq2, seq, cigar2.c_str(), 'G', 'A'));
                else
                    mismatch2 = EditDistance(refSeq2, seq, cigar2.c_str());
            }
        }
        PrintOutput(stdout, qname, aligned, cigLen, rname, pos, mismatch, gapOpen, gapExt, chr, startPos, mismatch2, gapOpen2, gapExt2, seq, cigar, refSeq, (char *) cigar2.c_str(), refSeq2);
    }
	fclose(stdin);
	fclose(stdout);
}


int main(int argc, char * argv[]) {
// Parsing Arugments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-r") == 0) {
            realPos = true;
            continue;
        }
        if (strcmp(argv[i], "-b") == 0) {
            bisSeq = true;
            continue;
        }
        if (strcmp(argv[i], "-s") == 0) {
            seqOutput = true;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0) {
            onlyPenalties = true;
            continue;
        }
        if (strcmp(argv[i], "-k") == 0) {
            keepSecondary = true;
            continue;
        }
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
        cerr << "Arguments:\n" <<
             "    -g <fasta file>              reference genome, mandatory\n" <<
             "    -i <SAM file>                the input SAM file to be analyzed, default: standard input\n" <<
             "    -o <result file>             the text file to which the results will be written, default: standard output\n" <<
             "    -s                           print SAM, CIGAR and alignment sequences in output\n" <<
             "    -b                           bisulfite-seq input, default: DNA-seq\n" <<
             "    -r                           use only if reads are generated via read_simul, or correct read location is provided in the read name in a similar format\n" <<
             "    -p                           only print the penalties (number of mismatchs, gap openings and gap extensions)\n" <<
             "    -k                           keep non-primary alignments and report them in the output\n";
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

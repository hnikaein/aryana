#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
using namespace std;

const long maxGenomeSize = 4e9;
const int maxChromosomeNum = 1000;
const int maxReadSize = 10000;
long chromPos[maxChromosomeNum];
unsigned long gs;
map <string, int> chromIndex;
map <int, string> chromName;
char * genome;
unsigned short * meth;
int chromNum = 0;
string genomeFile, annotationFile, outputFile, methInFile, methOutFile;
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
unsigned long long int n = 0, ni = 0, snp = 0, readl = 0;
double island = 0, cpg = 0, noncpg = 0, err = 0;
bool repeatMask = false, neg = false, pcr = false; // Mask repeat regions, produce reads from negative strand, produce PCR amplified reads

void ReadGenome(string genomeFile) {
    cerr << "Allocating memory..." << endl;
    genome = new char[maxGenomeSize];
    meth = new unsigned short[maxGenomeSize];
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
        int n = fscanf(f, "%s\n", fLine);
        if (n == EOF) break;
        n = strlen(fLine);
        if (fLine[0] == '>') {
            chromPos[chromNum++] = gs;
            if (chromNum > 1) cerr << " Length: " << chromPos[chromNum - 1] - chromPos[chromNum - 2] <<  endl;
            cerr << fLine;
            string name = fLine;
            name = name.substr(1, name.find("|")-1);
            chromIndex[name] = chromNum - 1;
            chromName[chromNum - 1] = name;
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
    for (long long i = 0; i < strlen(seq); i++)
        if (seq[i] == 'N' || seq[i] == 'n') {
            result = false;
            break;
        }
    return result;
}

bool ConvertPosition(int chr, long long wStart, long long wEnd, long long &start, long long &end)
{
    if (wStart < 1 || wStart + chromPos[chr] - 1 >= chromPos[chr + 1]) return false;

    start = wStart + chromPos[chr] - 1;
    end = wEnd + chromPos[chr] - 1;
    return true;
}

bool ConvertSequence(long chr, long long wStart, long long wEnd, double noncpg, double cpg) {
    //pos--; // Converting 1 base position to 0 base
    if (wStart < 1 || wStart + chromPos[chr] - 1 >= chromPos[chr + 1]) return false;

    //cout << chr << "\t" << wStart << "\t" << wEnd << endl;
    for (long long i = wStart + chromPos[chr] - 1; i < wEnd + chromPos[chr]; i++) {
        //cout << genome[i];
        if (genome[i] == 'c' || genome[i] == 'C') {
            if (genome[i+1] == 'g' || genome[i+1] == 'G') { // A CpG
                if ((double) rand () / RAND_MAX > cpg) // Non-methylated
                    genome[i] = 'T';
            } else {
                if ((double) rand() / RAND_MAX > noncpg)
                    genome[i] = 'T'; // Simulating bisulfite treatment for non-CpG cytosines
            }
        }
    }
    //cout << endl;
    return true;
}

void ConvertWholeGenome() {
    for (int i = 0; i < chromNum; i++)
        ConvertSequence(i, chromPos[i] + 1, chromPos[i+1], noncpg, cpg);
}

void AssignMethylationRatio(string cpgIslandFile) {
    unsigned short cpg1 = (unsigned short) (0.5 + cpg * 255), noncpg1 = (unsigned short) (0.5 + noncpg * 255), island1 = (unsigned short) (0.5 + island * 255);
    if (methInFile != "") {
        cerr << "Assigning methylation ratios based on file: " << methInFile << endl;
        ifstream f(methInFile.c_str());
        if (f.fail()) {
            cerr << "Error opening " << methInFile << endl;
            exit(1);
        }
        while (f.is_open() && ! f.eof()) {
            string chr;
            long long position;
            double m;
            f >> chr >> position >> m;
            position += chromPos[chromIndex[chr]] - 1;
            meth[position] = (unsigned short) (0.5 + m * 255);
        }
        f.close();
        return;
    }
    cerr << "Assignign methylation ratios" << endl;
    for (int i = 0; i < chromNum; i++)
        for (long long j = chromPos[i]; j < chromPos[i+1]; j++)
        {
            //cerr << i << "\t" << j << "\t" << genome[j] << "\t" << genome[j+1] << endl;
            if (genome[j] == 'c' || genome[j] == 'C')
                if ((j+1 < chromPos[i+1]) && (genome[j+1] == 'g' || genome[j+1] == 'G')) meth[j] = cpg1;
                else meth[j] = noncpg1;
            if (genome[j] == 'g' || genome[j] == 'G')
                if ((j > chromPos[i]) && (genome[j-1] == 'c' || genome[j-1] == 'C')) meth[j] = cpg1;
                else meth[j] = noncpg1;
        }


    if (cpgIslandFile == "") {
        cerr << "No CpG island location indicated. Proceeding without CpG island locations." << endl;
    } else {
        cerr << "Processing CpG island locations from file: " <<  annotationFile << endl;
        ifstream f(cpgIslandFile.c_str());
        if (! f.is_open()) {
            cerr << "Error: CpG island locations file not found or could not be opened" << endl;
            exit(1);
        }
        char fLine[10000], chrom[10], strand[10];
        long long wStart, wEnd, chr;
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
            chr = chromIndex[chrom];
            islands.push_back(window(chr, wStart, wEnd));
            for (long long i = wStart + chromPos[chr] - 1; i < wEnd + chromPos[chr]; i++) {
                if (genome[i] == 'c' || genome[i] == 'C')
                    if (i + 1 < gs && (genome[i+1] == 'g' || genome[i+1] == 'G')) meth[i] = island1;
                if (genome[i] == 'g' || genome[i] == 'G')
                    if (i > 0 && (genome[i-1] == 'c' || genome[i-1] == 'C')) meth[i] = island1;
            }
        }
    }
    //for (int i = 0; i < gs; i++) cerr << meth[i] << " ";
    if (methOutFile != "") {
        cerr << "Writing methylation ratio of the nucleotides in: " << methOutFile << endl;
        ofstream f(methOutFile.c_str());
        for (int i = 0; i < chromNum; i++)
            for (long long j = chromPos[i]; j < chromPos[i+1]; j++)
                if (meth[j] > 0) f << chromName[i] << "\t" << j+1-chromPos[i] << "\t" << setprecision(2) << (double)meth[j] / 255 << endl;
        f.close();
    }
}

void WriteGenome(string outputFile)
{
    FILE * f = fopen(outputFile.c_str(), "w");
    char s[100], tmp;
    for (int i = 0; i < chromNum; i++) {
        fprintf(f, "%s\n", chromName[i].c_str());
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

unsigned long long lrand() {
    return RAND_MAX * rand() + rand();
}

char Mutate(char base) {
    char b[4] = {'A','C','G','T'}, c;
    do {
        c = b[rand() % 4];
    } while (c == base || c+'a'-'A' == base);
    return c;
}

void revcomp(char * a, long long l) {
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
}

void PrintRead(FILE *f, int readNumber, int chr, long long p, char *quals) {
    char strand = '+', original='o';
    if (neg && (double) rand() / RAND_MAX < 0.5) strand='-';
    if (pcr && (double) rand() / RAND_MAX < 0.5) original = 'p';
    fprintf(f, "@%d|%s:%llu-%llu|%c%c\n", readNumber+1, chromName[chr].c_str(), p - chromPos[chr]+1, p-chromPos[chr] + readl, strand, original);
    char r[maxReadSize];
    memcpy(r, genome + p, readl);
    if (strand == '-') revcomp(r, readl);
    for (int i = 0; i < readl; i++) // Bisulfite conversion
        if ((r[i] == 'c' || r[i] == 'C') && (double) rand() * 255 / RAND_MAX >= meth[i+p]) r[i] = 'T';
    if (original == 'p') revcomp(r, readl);
    r[readl] = 0;
    for (int i = 0; i < readl; i++)
        if ((double) rand() / RAND_MAX < err) r[i] = Mutate(r[i]); // NGS Errors
    fprintf(f, "%-60s\n+\n%-60s\n", r,quals);
}

bool Repeat(const char *g, int length) {
    for (int i = 0; i < length; i++)
        if (g[i] >= 'a') return true;
    return false;
}

void SimulateReads(string outputFile) {
    FILE * f = fopen(outputFile.c_str(), "w");
    char quals[maxReadSize];
    memset(quals, 'I', readl);
    quals[readl] = 0;
    unsigned long long p;

    for (int i = 0; i < n; i++) {
        bool found;
        int j;
        do {
            found = true;
            p = lrand() % gs;
            for (j = 0; j <= chromNum; j++) // Checking the read to be inside a chromosome
                if (p < chromPos[j]) {
                    if (p+readl-1 >= chromPos[j]) found = false;
                    break;
                }
            if (found && memchr(genome + p, 'N', readl)) found = false;
            else if (repeatMask && Repeat(genome + p, readl)) found = false;
        } while (!found);
        j--;
        PrintRead(f, i, j, p, quals);
    }

    for (int i = 0; i < ni; i++) {
        bool found;
        int is;
        do {
            is = lrand() % islands.size();
            if (islands[is].wEnd - islands[is].wStart + 1 < readl) continue;
            p = rand() % (islands[is].wEnd - islands[is].wStart - readl + 2) + islands[is].wStart + chromPos[islands[is].chr];
            if (memchr(genome + p, 'N', readl)) continue;
            if (repeatMask && Repeat(genome + p, readl)) continue;
            break;
        } while (true);
        PrintRead(f, n+i, islands[is].chr, p, quals);
    }
    fclose(f);
}


void ProcessSNPs(void) {
    for (int i = 0; i < snp; i++) {
        unsigned long long j = lrand() % gs;
        genome[j] = Mutate(genome[j]);
    }
}



int main(int argc, char * argv[]) {
// Parsing Arugments
    if (argc < 6) {
        cerr << "Usage is -g <reference genome> -n <number of simulated reads from whole genome>  -ni <number of simulated reads from CpG-Islands>" << endl <<
             "-l <length of each read> -m (mask repeats) -mi <methylation ratio input file> -mo <methylation ratio output file> -ch <methylation ratio of CH dinucleotides>" << endl <<
             "-cg <methylation ratio of CG dinucleotides> -i <methylation ratio of CG dinucleotides in CpG-Islands> -e <sequencing errors ratio> -s <number of SNPs>"<<endl <<
             " -a <genomic map of CpG-Islands> -neg (generate reads from negative strand) -p (generate pcr amplified reads) -o <fastq output file>" << endl;
        exit(1);
    }
    for (int i = 1; i < argc; i++)
        if (i + 1 != argc) {
            if (strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
            else if (strcmp(argv[i], "-a") == 0)
                annotationFile = argv[++i];
            else if (strcmp(argv[i], "-o") == 0)
                outputFile = argv[++i];
            else if (strcmp(argv[i], "-n") == 0)
                n = atoi(argv[++i]);
            else if (strcmp(argv[i], "-ni") ==0)
                ni = atoi(argv[++i]);
            else if (strcmp(argv[i], "-ch")== 0)
                noncpg = atof(argv[++i]);
            else if (strcmp(argv[i], "-cg")==0)
                cpg = atof(argv[++i]);
            else if (strcmp(argv[i], "-i")==0)
                island = atof(argv[++i]);
            else if (strcmp(argv[i], "-e") == 0)
                err = atof(argv[++i]);
            else if (strcmp(argv[i], "-s") == 0)
                snp = atoi(argv[++i]);
            else if (strcmp(argv[i], "-l") == 0)
                readl = atoi(argv[++i]);
            else if (strcmp(argv[i], "-m") == 0)
                repeatMask = true;
            else if (strcmp(argv[i], "-neg") ==0)
                neg = true;
            else if (strcmp(argv[i], "-p") == 0)
                pcr = true;
            else if (strcmp(argv[i], "-mi") == 0)
                methInFile = argv[++i];
            else if (strcmp(argv[i], "-mo") == 0)
                methOutFile = argv[++i];
            else {
                cerr << "Not enough or invalid arguments."<< endl;
                exit(1);
            }
        }
    if (outputFile == "") {
        cerr << "No output file indicated. Casting output to BisSimul.fastq" << endl;
        outputFile = "BisSimul.fastq";
    }
// Reading input
    ReadGenome(genomeFile);
    ProcessSNPs();
    AssignMethylationRatio(annotationFile);
    SimulateReads(outputFile);
    delete [] genome;
}

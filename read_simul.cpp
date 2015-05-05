#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
using namespace std;

const int maxReadSize = 10000;
const int maxTries = 10000; // Maximum number of tries to find a suitable location for each read
vector<long long> chromPos, chromLen;
long long gs;
map <string, int> chromIndex;
map <int, string> chromName;
char * genome;
unsigned short * meth = 0;
unsigned int * totalCount = 0, *methylCount = 0;
int chromNum = 0;
string genomeFile, cpgIslandFile, methInFile, methOutFile, outputFileName, countOutFile;
FILE * outputFile = stdout, * outputFile2 = stdout;
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
long long int n = 1000, ni = 0, snp = 0, readl = 100;
double island = 0.1, cpg = 0.9, noncpg = 0.01, mismatchRate = 0, insRate = 0, delRate = 0;
bool bisSeq = false, repeatMask = false, neg = false, pcr = false, paired=false; // Mask repeat regions, produce reads from negative strand, produce PCR amplified reads, produce paired-end reads
int pairMinDis = 300, pairMaxDis = 1000; // Minimum and maximum distance between paired end reads
char bases[4] = {'A','C','G','T'};
void ReadGenome(string genomeFile) {
    cerr << "Allocating memory..." << endl;
    ifstream ifile(genomeFile.c_str());
    ifile.seekg(0, std::ios_base::end);	//seek to end
    //now get current position as length of file
    long long size = ifile.tellg();
    ifile.close();
    genome = new char[size];
    meth = new unsigned short[size];
    bzero(meth, size * sizeof(unsigned short));
    if (countOutFile != "") {
        totalCount = new unsigned int[size];
        methylCount = new unsigned int[size];
        bzero(totalCount, size * sizeof(unsigned int));
        bzero(methylCount, size * sizeof(unsigned int));
    }
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
//	for (int i = 0; i < chromNum; i++)
//		cerr << "ChromPos: " << chromPos[i] << " ChromLen: " << chromLen[i] << endl;
    cerr << " Length: " << chromPos[chromNum] - chromPos[chromNum - 1] <<  endl;
    fclose(f);
}

// Returns the chromosome number (zero based) if the given genome-wide location is completely inside a chromosome, or -1 otherwise

int InsideChromosome(long long start, long long end) {
    if (start < 0 || start > end) return -1;
    for (int j = 1; j <= chromNum; j++) // Checking the read to be inside a chromosome
        if (start < chromPos[j]) {
            if (end < chromPos[j]) return j - 1;
            else return -1; // Both start and end within the same chromosome
        }
    return -1; // Bigger than all chromosome locations
}

// Reads the CpG islands, checks their position to be inside chromosomes, and adds them to the islands vector

void ReadCpGIslands() {
    if (cpgIslandFile == "") {
        cerr << "No CpG island location indicated. Proceeding without CpG island locations." << endl;
        return;
    }
    cerr << "Processing CpG island locations from file: " <<  cpgIslandFile << endl;
    ifstream f(cpgIslandFile.c_str());
    if (! f.is_open()) {
        cerr << "Error: CpG island locations file not found or could not be opened" << endl;
        exit(1);
    }
    char fLine[10000], chrom[10];
    long long wStart, wEnd;
    while (f.good()) {
        fLine[0] = 0;
        f.getline(fLine, sizeof(fLine));
        if (! fLine[0]) continue;
        // cerr << fLine << endl;
        wStart = 0;
        sscanf(fLine, "%s %lld %lld", chrom, &wStart, &wEnd);
        if (! wStart || wStart >= wEnd || wEnd > chromLen[chromIndex[chrom]])
            cerr << "Omitting incorrect CpG island: " << chrom << ':' << wStart << '-' << wEnd << endl;
        else
            islands.push_back(window(chromIndex[chrom], wStart - 1, wEnd - 1));
    }
    f.close();
    cerr << "Read " << islands.size() << " CpG islands." << endl;
}

void AssignMethylationRatio(string cpgIslandFile) {
    unsigned short cpg1 = (unsigned short) (cpg * 100), noncpg1 = (unsigned short) (noncpg * 100), island1 = (unsigned short) (island * 100);
    char chr[1000];
    long long position;
    double m;
    if (methInFile != "") {
        cerr << "Assigning methylation ratios based on file: " << methInFile << endl;
        FILE * f = fopen(methInFile.c_str(), "r");
        if (!f) {
            cerr << "Error opening " << methInFile << endl;
            exit(1);
        }
        while (! feof(f)) {
            position=0;
            fscanf(f, "%s\t%lld\t%lf\n", chr, &position, &m);
            if (position <= 0) continue;
            position += chromPos[chromIndex[chr]] - 1;
            meth[position] = (unsigned short) (m * 100);
        }
        fclose(f);
        return;
    }
    cerr << "Assignign methylation ratios" << endl;
    for (int i = 0; i < chromNum; i++)
        for (long long j = chromPos[i]; j < chromPos[i+1]; j++)
        {
            //cerr << i << "\t" << j << "\t" << genome[j] << "\t" << genome[j+1] << endl;
            if (genome[j] == 'c' || genome[j] == 'C') {
                if ((j+1 < chromPos[i+1]) && (genome[j+1] == 'g' || genome[j+1] == 'G')) meth[j] = cpg1;
                else meth[j] = noncpg1;
            }
            if (genome[j] == 'g' || genome[j] == 'G') {
                if ((j > chromPos[i]) && (genome[j-1] == 'c' || genome[j-1] == 'C')) meth[j] = cpg1;
                else meth[j] = noncpg1;
            }
        }

    for (unsigned int j = 0; j < islands.size(); j++) {
        for (long long i = islands[j].wStart + chromPos[islands[j].chr]; i <= islands[j].wEnd + chromPos[islands[j].chr]; i++) {
            if (genome[i] == 'c' || genome[i] == 'C')
                if (i + 1 < gs && (genome[i+1] == 'g' || genome[i+1] == 'G')) meth[i] = island1;
            if (genome[i] == 'g' || genome[i] == 'G')
                if (i > 0 && (genome[i-1] == 'c' || genome[i-1] == 'C')) meth[i] = island1;
        }
    }
}

long long lrand() {
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
    delete [] b;
}

char RandBase() {
	return bases[rand() % 4];
}

// Just generates a single read, without header line but with quality line

void PrintSingleRead(FILE *f, long long p, char * quals, char strand, char original) {
    char r[maxReadSize], R[maxReadSize];
    memcpy(r, genome + p, readl);
    if (strand == '-') revcomp(r, readl);
	if (bisSeq) 
    for (int i = 0; i < readl; i++) // Bisulfite conversion
        if (r[i] == 'c' || r[i] == 'C') {
            if (strand == '+') {
                if (totalCount) totalCount[p + i]++;
                if ((double) rand() * 100.0 / RAND_MAX >= meth[i+p]) r[i] = 'T';
                else if (methylCount) methylCount[p+i]++;
            }
            else {
                if (totalCount) totalCount[p+readl-1-i]++;
                if ((double) rand() * 100.0 / RAND_MAX >= meth[p+readl-1-i]) r[i] = 'T';
                else  if (methylCount) methylCount[p+readl-1-i]++;
            }
        }
    if (original == 'p') revcomp(r, readl);
    r[readl] = 0;
	int RLen = 0;
    for (int i = 0; i < readl; i++) {
		double e = (double) rand() / RAND_MAX; 
        if (e < mismatchRate) R[RLen++] = Mutate(r[i]); // Mutation
		else if (e < insRate + mismatchRate) { // Insertion
			R[RLen++] = RandBase();
			i--;
		} else if (e < delRate + insRate + mismatchRate) { // Deletion, do nothing 
		} else R[RLen++] = r[i];
	}
	R[RLen] = 0;
	quals[RLen] = 0;
    fprintf(f, "%-60s\n+\n%-60s\n", R,quals);
	quals[RLen] = 'I';
}

// Can print both single and paired end reads, along with the header line

void PrintRead(int readNumber, int chr, long long p, char *quals, int pairDis) {
    char strand = '+', original='o';
    if (neg && (double) rand() / RAND_MAX < 0.5) strand='-';
    if (pcr && (double) rand() / RAND_MAX < 0.5) original = 'p';
    if (! paired) fprintf(outputFile, "@%d|%s:%llu-%llu|%c%c\n", readNumber+1, chromName[chr].c_str(), p - chromPos[chr]+1, p-chromPos[chr] + readl, strand, original);
    else {
        fprintf(outputFile, "@%d_1|%s:%llu-%llu|%c%c\n", readNumber+1, chromName[chr].c_str(), p - chromPos[chr]+1, p-chromPos[chr] + readl, strand, original);
        fprintf(outputFile2, "@%d_2|%s:%llu-%llu|%c%c\n", readNumber+1, chromName[chr].c_str(), p - chromPos[chr]+readl+pairDis+1, p-chromPos[chr] + 2 * readl + pairDis, strand, original);
    }
    PrintSingleRead(outputFile, p, quals, strand, original);
    if (paired) PrintSingleRead(outputFile2, p + readl + pairDis, quals, strand, original);
}

// Looks for a genomic region to find out any repeat regions (lower-case nucleotides)

bool Repeat(const char *g, int length) {
    for (int i = 0; i < length; i++)
        if (g[i] >= 'a') return true;
    return false;
}

// Returns the chromosome number if the given location is a good posotion for generating a read, or -1 othersise

int CheckPosition(long long p, int pairDis) {
    int chr = -1;
    if ((chr = InsideChromosome(p, p + readl - 1)) == -1) return -1;
    if (memchr(genome + p, 'N', readl)) return -1;
    if (repeatMask && Repeat(genome + p, readl)) return -1;
    if (paired) {
        if (InsideChromosome(p, p + 2 * readl + pairDis - 1) != chr) return -1;
        if (memchr(genome + p + readl + pairDis, 'N', readl)) return -1;
        if (repeatMask && Repeat(genome + p + readl + pairDis, readl)) return -1;
    }
    return chr;
}

// Simulates the reads

void SimulateReads() {
    char quals[maxReadSize];
    memset(quals, 'I', sizeof(quals));
    long long p;
    int pairDis = 0, chr;
    if (! gs) {
        cerr << "Error: the length of genome is zero." << endl;
        exit(-1);
    }

// Simulating normal reads
    for (int i = 0; i < n; i++) {
        int tries = 0;
        do {
            if (tries++ > maxTries) {
                cerr << "Error: Could not find any suitable location for reads after " << maxTries << " tries." << endl;
                exit(-1);
            }
            p = lrand() % gs;
            pairDis = lrand() % (pairMaxDis - pairMinDis + 1) + pairMinDis;
            if ((chr = CheckPosition(p, pairDis)) != -1) break;
        } while (true);
        PrintRead(i, chr, p, quals, pairDis);
    }

// Simulating CpG-island reads

    if (ni > 0 && islands.size() == 0) {
        cerr << "Error: There are no CpG islands read, while there should be reads produced from CpG islands." << endl;
        exit(-1);
    }
    for (int i = 0; i < ni; i++) {
        int tries = 0;
        int is;
        do {
            if (tries++ > maxTries) {
                cerr << "Error: Could not find any suitable location for CpG-island reads after " << maxTries << " tries." << endl;
                exit(-1);
            }
            is = lrand() % islands.size();
            p = rand() % (islands[is].wEnd - (islands[is].wStart - readl + 1) + 1) + islands[is].wStart - readl + 1 + chromPos[islands[is].chr]; // To have at least 1bp overlap between read and CpG-island
//			cerr << islands.size() << " " << is << " " << p << " " << islands[is].chr << " " << CheckPosition(p, pairDis) << endl;
            if ((chr = CheckPosition(p, pairDis)) == islands[is].chr) break;
        } while (true);
        PrintRead(n+i, chr, p, quals, pairDis);
    }

// Closing output
    if (outputFile != stdout) fclose(outputFile);
    if (paired && outputFile2 != stdout) fclose(outputFile2);
}


void ProcessSNPs(void) {
    if (! gs) {
        cerr << "Error: the length of genome is zero." << endl;
        exit(-1);
    }
    for (int i = 0; i < snp; i++) {
        long long j = lrand() % gs;
        genome[j] = Mutate(genome[j]);
    }
}

void WriteMethylationRatios() {
    if (methOutFile != "") {
        cerr << "Writing methylation ratio of the nucleotides in: " << methOutFile << endl << "Please be patient, this might take some time depending on the genome size" << endl;
        FILE * f = fopen(methOutFile.c_str(), "w");
        fprintf(f, "chr\tpos\tmeth\n");
        for (int i = 0; i < chromNum; i++)
            for (long long j = chromPos[i]; j < chromPos[i+1]; j++)
                if (meth[j] > 0) fprintf(f, "%s\t%lld\t%.2f\n", chromName[i].c_str(), j+1-chromPos[i], (double) meth[j] / 100.0);
        fclose(f);
    }
}

void WriteReadCounts() {
    if (countOutFile != "") {
        cerr << "Writing the number of reads covering each cytosine (converted, unconverted, ratio) in: " << countOutFile << endl << "Please be patient, this might take some time depending on the genome size" << endl;
        FILE * f = fopen(countOutFile.c_str(), "w");
        fprintf(f, "chr\tpos\tnTotal\tnMeth\tmeth\n");
        for (int i = 0; i < chromNum; i++)
            for (long long j = chromPos[i]; j < chromPos[i+1]; j++)
                if (totalCount[j] > 0) fprintf(f, "%s\t%lld\t%d\t%d\t%.2f\n", chromName[i].c_str(), j+1-chromPos[i], totalCount[j], methylCount[j], (double) methylCount[j] / totalCount[j]);
        fclose(f);
    }
}

int main(int argc, char * argv[]) {
// Parsing Arugments
    for (int i = 1; i < argc; i++) {
        // The single arguments
        if (strcmp(argv[i], "-m") == 0) {
            repeatMask = true;
            continue;
        }
        if (strcmp(argv[i], "-neg") ==0) {
            neg = true;
            continue;
        }
        if (strcmp(argv[i], "-p") == 0) {
            pcr = true;
            continue;
        }
        if (strcmp(argv[i], "-P") == 0) {
            paired = true;
            continue;
        }

        if (i < argc - 1) { // The paired arguments
            if (strcmp(argv[i], "-g") == 0)
                genomeFile = argv[++i];
            else if (strcmp(argv[i], "-a") == 0)
                cpgIslandFile = argv[++i];
            else if (strcmp(argv[i], "-o") == 0)
                outputFileName = argv[++i];
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
            else if (strcmp(argv[i], "-rm") == 0)
                mismatchRate = atof(argv[++i]);
			else if (strcmp(argv[i], "-ri") == 0)
				insRate = atof(argv[++i]);
			else if (strcmp(argv[i], "-rd") == 0)
				delRate = atof(argv[++i]);
			else if (strcmp(argv[i], "-b") == 0)
				bisSeq = true;
            else if (strcmp(argv[i], "-s") == 0)
                snp = atoi(argv[++i]);
            else if (strcmp(argv[i], "-l") == 0)
                readl = atoi(argv[++i]);
            else if (strcmp(argv[i], "-mi") == 0)
                methInFile = argv[++i];
            else if (strcmp(argv[i], "-mo") == 0)
                methOutFile = argv[++i];
            else  if (strcmp(argv[i], "-co") == 0)
                countOutFile = argv[++i];
            else if (strcmp(argv[i], "-Pmin") == 0)
                pairMinDis = atoi(argv[++i]);
            else if (strcmp(argv[i], "-Pmax") == 0)
                pairMaxDis = atoi(argv[++i]);
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
        cerr << "Arguments:" << endl <<
				"-g     <reference genome, mandatory>" << endl <<
				"-n     <number of simulated reads from whole genome, default=1000>" << endl <<
			    "-l     <length of each read, default=100>" << endl << 
				"-o     <fastq output file without extension, default=stdout>" << endl <<
				"-m     (mask repeats)" << endl <<
				"-rm    <sequencing mismatch rate, between 0 and 1, default=0>" << endl << 
				"-ri    <sequencing insertion rate, between 0 and 1, default=0>" << endl <<
				"-rd    <sequencing deletion rate, between 0 and 1, default=0>" << endl <<
				"-s     <number of SNPs, default=0>" << endl <<		
			    "-neg   (simulate reads from negative strand)" << endl <<  
				"-p     (simulate PCR amplified reads)" << endl << endl <<
				"Paired end arguments:" << endl <<
				"-P     (produce paired-end reads)" << endl <<
				"-Pmin  <min distance between a pair of reads,default=300>" << endl <<
				"-Pmax  <max distance between a pair of reads, default=1000>" << endl << endl <<
			    "Bisulfite-Sequencing arguments:" << endl <<
				"-b     (produce bis-seq reads)" << endl <<
				"-ni    <number of simulated reads from CpG-Islands, default=0>" << endl <<
             	"-a     <genomic map of CpG-Islands, mandatory if -ni greater than 0>" << endl << 
             	"-mi    <methylation ratio input file>" << endl <<
				"-mo    <methylation ratio output file>" << endl << 
				"-co    <count of reads covering each cytosine output file>" << endl <<
             	"-ch    <methylation ratio of CH dinucleotides, default=0.01>" << endl << 
				"-cg    <methylation ratio of CG dinucleotides outside CpG-Islands, default=0.9>" << endl <<
             	"-i     <methylation ratio of CG dinucleotides in CpG-Islands, default=0.1>" << endl;
        exit(1);
    }

    if (pairMinDis > pairMaxDis || pairMinDis < 0) {
        cerr << "Error: Min distance between pair of reads should be a non-negative integer, smaller than max distance between them" << endl;
        exit(1);
    }
	if (insRate + delRate + mismatchRate > 1) {
		cerr << "Error: Sum of insertion, deletion and mismatch ratios should be between 0 and 1" << endl;
		exit(1);
	}
// Reading input
    ReadGenome(genomeFile);
    ReadCpGIslands();
    ProcessSNPs();
	if (bisSeq)
	    AssignMethylationRatio(cpgIslandFile);
    if (outputFileName != "") {
        if (paired) {
            outputFile = fopen((outputFileName + "_1.fastq").c_str(), "w");
            outputFile2 = fopen((outputFileName + "_2.fastq").c_str(), "w");
        } else outputFile = fopen((outputFileName + ".fastq").c_str(), "w");
    }
    cerr << "Simulating reads..." << endl;
    SimulateReads();
	if (bisSeq)
	    WriteMethylationRatios();
    WriteReadCounts();
    cerr << "Finished." << endl;
    delete [] genome;
    delete [] meth;
    if (totalCount) {
        delete[] totalCount;
        delete[] methylCount;
    }
}

#define __STDC_FORMAT_MACROS

#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cmath>
#include <cctype>
#include "bwt.h"
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <iostream>

#include "const.h"

using namespace std;

const int unalignedPenalty = 1000000;
const int discordPenalty = 100000;
unsigned long long gs;
vector<string> chromName;
char *genome;
int chromNum, minD = 0, maxD = 10000;
FILE *outputFile = nullptr, *headerFile = nullptr;
bool paired = false;
int chosen[BS_GENOMES_COUNT];

enum orientation_t {
    fr, rf, ff, all
} orientation = all;

struct Island {
    int chr;
    long long start, end;

    Island(int c, long long s, long long e) {
        chr = c;
        start = min(s, e);
        end = max(s, e);
    }

    bool operator<(const Island &val) const {
        return chr < val.chr || (chr == val.chr && end < val.start);
    }
};

vector<Island> islands;
int highPenalty = 0, lowPenalty = 0;
long readPenalties[BS_GENOMES_COUNT];

char *samNames[BS_GENOMES_COUNT];
FILE *samFiles[BS_GENOMES_COUNT];
char str[MAX_READ_LENGTH];


struct Chrom {
    string chrName;
    long long chrStart;
}; // for converting refrence position to chromosome based position

vector<Chrom> chrom;


// Reads genome and saves start position of chromosomes
int ReadGenome(char *genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info;
    FILE *f;
    if (stat(genomeFile, &file_info) == -1) {
        fprintf(stderr, "Could not get the information of file %s\nplease make sure the file exists\n", genomeFile);
        return -1;
    }

    f = fopen(genomeFile, "r");
    if (!f) {
        fprintf(stderr, "Error opening reference file: %s\n", genomeFile);
        exit(-1);
    }
    off_t file_size_bytes = file_info.st_size;
    long long reference_size = 1 + ceil(((double) file_size_bytes) / (double) (sizeof(char)));
    genome = (char *) malloc(reference_size * sizeof(char));
    memset(genome, 0, reference_size * sizeof(char));
    gs = 0;
    chromNum = 0;
    fprintf(stderr, "Reading genome...\n");
    char fLineMain[1000000];
    Chrom ch;
    while (!feof(f)) {
        if (!fgets(fLineMain, sizeof(fLineMain), f)) break;
        int n = strlen(fLineMain), start = 0;
        while (n > 0 && fLineMain[n - 1] <= ' ') n--;
        fLineMain[n] = 0;
        while (start < n && fLineMain[start] <= ' ') start++;
        if (start >= n) continue;
        char *fLine = fLineMain + start;
        n -= start;

        if (fLine[0] == '>') {
            ch.chrStart = gs;
            string name = fLine;
            if (name.find(" ") != string::npos) name = name.substr(1, name.find(" ") - 1);
            else name = name.substr(1, name.size() - 1);
            //cerr << name << endl;
            ch.chrName = name;
            chrom.push_back(ch);
            chromNum++;
        } else {
            memcpy(genome + gs, fLine, n);
            gs += n;
        }
    }
    fclose(f);
    return 0;
}


// returns minimum penalty from scores for each genome , if penalties are equal , returns original genome
int min_penalty() {
    int j = 0;
    long min = readPenalties[0];
    for (int i = 1; i < BS_GENOMES_COUNT; i++) {
        // cerr<<readPenalties[i]<<endl;
        if (min > readPenalties[i]) {
            min = readPenalties[i];
            j = i;
        }
    }
    return j;
}

char getNuc(uint64_t place) {
    return genome[place];
}

inline void ToLower(char *s) {
    for (unsigned int i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z')
            s[i] += 'a' - 'A';
}


int ChromIndex(char *chr) {
    for (int i = 0; i < chrom.size(); i++) {
        if (chrom[i].chrName == chr)
            return i;
    }
    return -1;

}


void ReadCpGIslands(char *annotationFile) {
    FILE *file = fopen(annotationFile, "r"); /* should check the result */
    char fLine[10000], chrom_name[MAX_CHR_NAME_LENGTH];
    unsigned long long start, end;
    do {
        if (nullptr == fgets(fLine, sizeof(fLine), file)) {
            fclose(file);
            break;
        }
        sscanf(fLine, "%s %llu %llu", chrom_name, &start, &end);
        islands.emplace_back(ChromIndex(chrom_name), start, end);
    } while (true);
    sort(islands.begin(), islands.end());
}

bool isInIsland(int chr, long long pos) {
    return binary_search(islands.begin(), islands.end(), Island(chr, pos, pos));
}


void CalcPenalties(uint64_t ref_i, char read, long readNum, int chr, uint64_t chrPos, int flag, int flag2) {
    read = toupper(read);
    genome[ref_i] = toupper(genome[ref_i]);
    genome[ref_i + 1] = toupper(genome[ref_i + 1]);
    genome[ref_i - 1] = toupper(genome[ref_i - 1]);

    if (flag2) {
        if (read == 'A' || read == 'G') {
            if (getNuc(ref_i) != read)
                readPenalties[readNum] += highPenalty;
        } else { //read = C or T

            if (read == 'T' && getNuc(ref_i) == 'C') {
                if (getNuc(ref_i + 1) == 'G') { // in the CpG context
                    if (isInIsland(chr, chrPos)) // in CpG and also island
                        readPenalties[readNum] += 0;
                    else
                        readPenalties[readNum] += highPenalty;

                } else // out of CpG context

                    readPenalties[readNum] += lowPenalty;

            } else if (read == 'C' && getNuc(ref_i) == 'C') {
                if (getNuc(ref_i + 1) == 'G') { // in the CpG context
                    int temp = isInIsland(chr, chrPos);
                    if (temp == 1)  // in CpG and also island
                        readPenalties[readNum] += highPenalty;
                    else            // in CpG and out of island
                        readPenalties[readNum] += lowPenalty;

                } else              // out of CpG context
                    readPenalties[readNum] += highPenalty;

            } else if (read == 'C' && getNuc(ref_i) == 'T')
                readPenalties[readNum] += highPenalty;

        }
    } else {
        if (read == 'C' || read == 'T') {
            if (getNuc(ref_i) != read)
                readPenalties[readNum] += highPenalty;
        } else { //read = G or A

            if (read == 'A' && getNuc(ref_i) == 'G') {
                if (getNuc(ref_i - 1) == 'C') { // in the CpG context
                    if (isInIsland(chr, chrPos)) { // in CpG and also island
                        readPenalties[readNum] += 0;
                    } else {
                        readPenalties[readNum] += highPenalty;
                    }
                } else { // out of CpG context

                    readPenalties[readNum] += lowPenalty;
                }
            } else if (read == 'G' && getNuc(ref_i) == 'G') {
                if (getNuc(ref_i - 1) == 'C') { // in the CpG context
                    int temp = isInIsland(chr, chrPos);
                    if (temp == 1) { // in CpG and also island
                        readPenalties[readNum] += highPenalty;
                    } else {
                        readPenalties[readNum] += lowPenalty;
                    }
                } else { // out of CpG context
                    readPenalties[readNum] += highPenalty;
                }

            } else if (read == 'G' && getNuc(ref_i) == 'A')
                readPenalties[readNum] += highPenalty;
        }

    }

}

// Reads Cigar sequence and call CalcPenalties for calculating penalties for each base
void readCigar(char *cigar, uint64_t ref_i, char *seq_string, long long readNum, int chr, uint64_t chrPos, int flag,
               int refNum) {
    int pos = 0;
    int value = 0;
    uint64_t ref_index = ref_i;
    long read_index = 0;
    strcpy(str, seq_string);
    while (true) {
        if (!isdigit(cigar[pos])) {
            if (value > 0) {
                if (cigar[pos] == 'M' || cigar[pos] == 'm') {
                    int j;
                    for (j = 0; j < value; j++) {
                        CalcPenalties(ref_index, seq_string[read_index], readNum, chr, chrPos, flag, refNum);
                        ref_index++;
                        read_index++;
                    }
                } else if (cigar[pos] == 'D' || cigar[pos] == 'd') {
                    ref_index += value;
                    readPenalties[readNum] += highPenalty * value;      //high penalty for insertion
                } else if (cigar[pos] == 'I' || cigar[pos] == 'i') {
                    read_index += value;
                    readPenalties[readNum] += highPenalty * value;      //high penalty for deletion
                }
            } else if (cigar[pos] == '*') {
                readPenalties[readNum] += unalignedPenalty;     // maximum penalty for not aligned reads
                break;
            }
            if (cigar[pos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[pos] - '0';
        }
        pos++;
    }
}

void WriteHeader() {
    // Printing header of original samfile to standard output
    char line[MAX_SAM_LINE_LENGTH];
    string header;
    FILE *samFile = fopen(samNames[0], "r");
    if (fgets(line, 10000, samFile))
        while (line[0] == '@') {
            header += line;
            if (!fgets(line, 10000, samFile))
                break;
        }
    fprintf(headerFile, "%s", header.c_str());
    fclose(samFile);
    if (headerFile != outputFile) fclose(headerFile);
}

bool ReadLine(char *line, int len, FILE *f) {
    do {
        if (fgets(line, len, f) == nullptr) return false;
    } while (line[0] == '@');
    return true;
}

void Process() {
    for (int i = 0; i < BS_GENOMES_COUNT; i++) samFiles[i] = fopen(samNames[i], "r");
    char *rname[BS_GENOMES_COUNT], *cigar[BS_GENOMES_COUNT], *rname2[BS_GENOMES_COUNT], *cigar2[BS_GENOMES_COUNT];
    int flag[BS_GENOMES_COUNT], flagTwo[BS_GENOMES_COUNT];
    uint64_t pos[BS_GENOMES_COUNT], pos2[BS_GENOMES_COUNT];
    uint32_t mapq[BS_GENOMES_COUNT], mapq2[BS_GENOMES_COUNT];
    long long int tlen, tlen2;
    for (int i = 0; i < BS_GENOMES_COUNT; i++) {
        rname[i] = new char[MAX_CHR_NAME_LENGTH];
        cigar[i] = new char[MAX_READ_LENGTH * 2];
        rname2[i] = new char[MAX_CHR_NAME_LENGTH];
        cigar2[i] = new char[MAX_READ_LENGTH * 2];
    }
    char qname[10000], rnext[10000], pnext[10000], seq_string[MAX_READ_LENGTH], quality_string[MAX_READ_LENGTH]; // (Ali) please double check the limitations for qname, rnext, pnext and copy
    char qname2[10000], rnext2[10000], pnext2[10000], seq_string2[MAX_READ_LENGTH], quality_string2[MAX_READ_LENGTH];
    char line[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH], line2[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH];
    bzero(chosen, sizeof(chosen));

    while (true) {
        for (int i = 0; i < BS_GENOMES_COUNT; i++) {
            readPenalties[i] = 0;
            if (!ReadLine(line[i], MAX_SAM_LINE_LENGTH, samFiles[i])) return;
            sscanf(line[i], "%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", qname, &flag[i], rname[i],
                   &pos[i], &mapq[i], cigar[i], rnext, pnext, &tlen, seq_string, quality_string);
            int index = ChromIndex(rname[i]);
            if (index == -1) readPenalties[i] += unalignedPenalty;
            else {
                int flag2 = 0;
                if (i <= 2) // if read has been aligned to first four genomes
                    flag2 = 1;
                if (i > 2) // if read has been aligned to last three genome
                    flag2 = 0;
                readCigar(cigar[i], pos[i] + chrom[index].chrStart - 1, seq_string, i, index, pos[i], flag[i], flag2);
            }

            if (paired) {
                if (!ReadLine(line2[i], MAX_SAM_LINE_LENGTH, samFiles[i])) {
                    cerr << "Error: unexpected end of input SAM file while reading the paired end.\n";
                    return;
                }
                sscanf(line2[i], "%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", qname2, &flagTwo[i],
                       rname2[i], &pos2[i], &mapq2[i], cigar2[i], rnext2, pnext2, &tlen2, seq_string2, quality_string2);
                int index2 = ChromIndex(rname2[i]);
                if (index2 == -1) readPenalties[i] += unalignedPenalty;
                else {
                    int flag2 = 0;
                    if (i <= 2) // if read has been aligned to first four genomes
                        flag2 = 1;
                    if (i > 2) // if read has been aligned to last three genome
                        flag2 = 0;
                    readCigar(cigar2[i], pos2[i] + chrom[index2].chrStart - 1, seq_string2, i, index2, pos2[i],
                              flagTwo[i], flag2);
                }
                orientation_t orien;
                if (index >= 0 && index2 >= 0) {
                    if ((flag[i] & 16) == (flagTwo[i] & 16)) orien = ff;
                    else if (((flag[i] & 16) == 0 && pos[i] < pos2[i]) ||
                             ((flag[i] & 16) != 0 && pos[i] > pos2[i]))
                        orien = fr;
                    else orien = rf;
                    if (index != index2 || llabs((signed) (pos[i] - pos2[i])) < minD ||
                        llabs((signed) (pos[i] - pos2[i])) > maxD || (orientation != all && orientation != orien))
                        readPenalties[i] += discordPenalty;
                }
            }
            cerr << i << "\t" << readPenalties[i] << endl;
        }
        int min = min_penalty();
        chosen[min]++; // shows how many times a genome has been selected
        fprintf(outputFile, "%s", line[min]);
        if (paired) fprintf(outputFile, "%s", line2[min]);
    }
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr,
                "Usage: align_bs -x <reference genome> -s <input sam files> -c <CpG-island map> -p <penalties> [-o <output SAM file>] [-h <SAM header file>] [-l/--limit <maximum allowed mismatches>]\n");
        fprintf(stderr, "Additional arguments for paired-end reads:\n");
        fprintf(stderr,
                "    -P to consider the input as paired-end, -m/--min for minimum and -M/--max for maximum distance, and --fr/rf/ff for orientation.\n");
        return -1;
    }
    static struct option long_options[] = {
            {"samfiles",    required_argument, nullptr, 's'},
            {"reference",   required_argument, nullptr, 'x'},
            {"CpG-islands", required_argument, nullptr, 'c'},
            {"penalties",   required_argument, nullptr, 'p'},
            {"output",      required_argument, nullptr, 'o'},
            {"header",      required_argument, nullptr, 'h'},
            {"paired-end",  no_argument,       nullptr, 'U'},
            {"min",         required_argument, nullptr, 'm'},
            {"max",         required_argument, nullptr, 'M'},
            {"fr",          no_argument,       nullptr, 1},
            {"rf",          no_argument,       nullptr, 2},
            {"ff",          no_argument,       nullptr, 3},
            {"limit",       required_argument, nullptr, 'l'}
    };
    int option_index = 0;
    char *samFileName = 0, *tmp, c;
    char *annotationFile, *referenceName;
    outputFile = stdout;
    while ((c = getopt_long(argc, argv, "x:s:c:p:o:h:Pm:M:123l:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'x':
                referenceName = strdup(optarg);
                break;
            case 's':
                samFileName = strdup(optarg);
                tmp = strdup(samFileName);
                for (int i = 0; i < BS_GENOMES_COUNT; i++) {
                    samNames[i] = (char *) malloc(strlen(samFileName) + 10);
                    sprintf(samNames[i], "%s-%d", samFileName, i);
                }
                break;
            case 'c':
                annotationFile = (char *) malloc(strlen(optarg) + 1);
                strcpy(annotationFile, optarg);
                break;
            case 'p':
                lowPenalty = atoi(optarg);
//            medPenalty = atoi(argv[optind]);
                highPenalty = atoi(argv[optind + 1]);
                optind = optind + 2;
                break;
            case 'o':
                outputFile = fopen(optarg, "w");
                break;
            case 'h':
                headerFile = fopen(optarg, "w");
                break;
            case 'P':
                paired = true;
                break;
            case 'm':
                minD = atoi(optarg);
                break;
            case 'M':
                maxD = atoi(optarg);
                break;
            case 1:
                orientation = fr;
                break;
            case 2:
                orientation = rf;
                break;
            case 3:
                orientation = ff;
                break;
            case 'l':
                // limit = atoi(optarg);
                break;
        }
    }
    if (!headerFile) headerFile = outputFile;

    if (!samFileName) {
        cerr << "Error: input SAM file name not found." << endl;
        return 1;
    }
    ReadGenome(referenceName);
    ReadCpGIslands(annotationFile);
    WriteHeader();
    Process();
    for (int j = 0; j < BS_GENOMES_COUNT; j++)
        fprintf(stderr, "Number of reads aligned to genome %d: %d\n", j, chosen[j]);

    if (outputFile != stdout) fclose(outputFile);
}


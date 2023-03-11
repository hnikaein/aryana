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

#define UNALIGNED_PENALTY  1000000
#define DISCORD_PENALTY  100000

unsigned long long gs;
vector<string> chromName;
char *genome;
int chromNum, paired_min_distance = 0, paired_max_distance = 10000;
FILE *outputFile = nullptr, *headerFile = nullptr;
bool paired = false;
int chosen[BS_GENOMES_COUNT];
char *strtox_temp;

enum orientation_t {
    fr, rf, ff, all
} orientation = all;

struct Island {
    int chr;
    uint64_t start, end;

    Island(int c, uint64_t s, uint64_t e) {
        chr = c;
        start = min(s, e);
        end = max(s, e);
    }

    bool operator<(const Island &val) const {
        return chr < val.chr || (chr == val.chr && end < val.start);
    }
};

vector<Island> islands;
int mismatch_penalty = 0, go_penatly = 0, ge_penalty = 0;

char *samNames[BS_GENOMES_COUNT];
FILE *samFiles[BS_GENOMES_COUNT];
char str[MAX_READ_LENGTH];

int readPenalties[2 * BS_GENOMES_COUNT];
char rname[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH], rname2[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH];
int flag[BS_GENOMES_COUNT], flagTwo[BS_GENOMES_COUNT];
uint64_t pos[BS_GENOMES_COUNT], pos2[BS_GENOMES_COUNT];


struct Chrom {
    string chrName;
    unsigned long long chrStart{};
}; // for converting refrence position to chromosome based position

vector<Chrom> chrom;


// Reads genome and saves start position of chromosomes
int ReadGenome(char *genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info{};
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
    auto reference_size = (long long) (1 + ceil(((double) file_size_bytes) / (double) (sizeof(char))));
    genome = (char *) malloc(reference_size * sizeof(char));
    memset(genome, 0, reference_size * sizeof(char));
    gs = 0;
    chromNum = 0;
    fprintf(stderr, "Reading genome...\n");
    char fLineMain[1000000];
    Chrom ch;
    while (!feof(f)) {
        if (!fgets(fLineMain, sizeof(fLineMain), f)) break;
        int n = (int) strlen(fLineMain), start = 0;
        while (n > 0 && fLineMain[n - 1] <= ' ') n--;
        fLineMain[n] = 0;
        while (start < n && fLineMain[start] <= ' ') start++;
        if (start >= n) continue;
        char *fLine = fLineMain + start;
        n -= start;

        if (fLine[0] == '>') {
            ch.chrStart = gs;
            string name = fLine;
            if (name.find(' ') != string::npos) name = name.substr(1, name.find(' ') - 1);
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


inline char getNuc(uint64_t place) {
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

void find_min_penalties(int *min, int *min2) {
    *min = *min2 = 0;
    int min_penalty = readPenalties[0] + (paired ? readPenalties[0 + BS_GENOMES_COUNT] : 0);
    for (int i = 0; i < BS_GENOMES_COUNT; i++) {
        if (!paired) {
            if (min_penalty > readPenalties[i]) {
                min_penalty = readPenalties[i];
                *min = i;
            }
        } else {
            for (int j = 0; j < BS_GENOMES_COUNT; j++) {
                int new_penalty = readPenalties[i] + readPenalties[j + BS_GENOMES_COUNT];
                orientation_t orien;
                int index = ChromIndex(rname[i]);
                int index2 = ChromIndex(rname2[i]);
                if (index >= 0 && index2 >= 0) {
                    if ((flag[i] & 16) == (flagTwo[j] & 16))
                        orien = ff;
                    else if (((flag[i] & 16) == 0 && pos[i] < pos2[j]) ||
                             ((flag[i] & 16) != 0 && pos[i] > pos2[j]))
                        orien = fr;
                    else
                        orien = rf;
                    if (index != index2 || llabs((signed) (pos[i] - pos2[j])) < paired_min_distance ||
                        llabs((signed) (pos[i] - pos2[j])) > paired_max_distance ||
                        (orientation != all && orientation != orien))
                        new_penalty += DISCORD_PENALTY;
                }
                if (min_penalty > new_penalty) {
                    min_penalty = new_penalty;
                    *min = i;
                    *min2 = j;
                }
            }
        }
    }
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

bool isInIsland(int chr, uint64_t chr_pos) {
    return binary_search(islands.begin(), islands.end(), Island(chr, chr_pos, chr_pos));
}


void CalcPenalties(uint64_t ref_i, char read, long read_penalties_id, int chr, uint64_t chrPos, int refCT) {
    read = (char) toupper(read);
    genome[ref_i] = (char) toupper(genome[ref_i]);
    genome[ref_i + 1] = (char) toupper(genome[ref_i + 1]);
    genome[ref_i - 1] = (char) toupper(genome[ref_i - 1]);
    char methyl_nuc = 'G', bis_methyl_nuc = 'A', CPG_next_nuc = 'C';
    int next_nuc = -1;
    if (refCT) {
        methyl_nuc = 'C';
        bis_methyl_nuc = 'T';
        next_nuc = 1;
        CPG_next_nuc = 'G';
    }

    if ((read != bis_methyl_nuc && read != methyl_nuc) || (getNuc(ref_i) != methyl_nuc)) {
        if (getNuc(ref_i) != read)
            readPenalties[read_penalties_id] += mismatch_penalty;
    } else {
        bool is_methylated = read != bis_methyl_nuc;
        if (getNuc(ref_i + next_nuc) == CPG_next_nuc && !isInIsland(chr, chrPos)) // in CpG and out of island
            readPenalties[read_penalties_id] += is_methylated ? 0 : mismatch_penalty;
        else
            readPenalties[read_penalties_id] += is_methylated ? mismatch_penalty : 0;
    }
}

// Reads Cigar sequence and call CalcPenalties for calculating penalties for each base
void readCigar(char *cigar, uint64_t ref_i, char *seq_string, int read_penalties_id, int chr, uint64_t chrPos,
               int refCT) {
    int cpos = 0;
    int value = 0;
    uint64_t ref_index = ref_i;
    long read_index = 0;
    strcpy(str, seq_string);
    while (true) {
        if (!isdigit(cigar[cpos])) {
            if (value > 0) {
                if (cigar[cpos] == 'M' || cigar[cpos] == 'm') {
                    int j;
                    for (j = 0; j < value; j++) {
                        CalcPenalties(ref_index, seq_string[read_index], read_penalties_id, chr, chrPos, refCT);
                        ref_index++;
                        read_index++;
                    }
                } else if (cigar[cpos] == 'D' || cigar[cpos] == 'd') {
                    ref_index += value;
                    readPenalties[read_penalties_id] += ge_penalty * (value - 1) + go_penatly;      //high penalty for insertion
                } else if (cigar[cpos] == 'I' || cigar[cpos] == 'i') {
                    read_index += value;
                    readPenalties[read_penalties_id] += ge_penalty * (value - 1) + go_penatly;      //high penalty for deletion
                }
            } else if (cigar[cpos] == '*') {
                readPenalties[read_penalties_id] += UNALIGNED_PENALTY;     // maximum penalty for not aligned reads
                break;
            }
            if (cigar[cpos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[cpos] - '0';
        }
        cpos++;
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

void Process(bool print_output) {
    for (int i = 0; i < BS_GENOMES_COUNT; i++) samFiles[i] = fopen(samNames[i], "r");
    char line[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH], cigar[BS_GENOMES_COUNT][MAX_READ_LENGTH * 2],
            qname[MAX_READ_NAME_LENGTH], rnext[MAX_CHR_NAME_LENGTH], pnext[30], seq_string[MAX_READ_LENGTH], quality_string[MAX_READ_LENGTH],
            line2[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH], cigar2[BS_GENOMES_COUNT][MAX_READ_LENGTH * 2],
            qname2[MAX_READ_NAME_LENGTH], rnext2[MAX_CHR_NAME_LENGTH], pnext2[30], seq_string2[MAX_READ_LENGTH], quality_string2[MAX_READ_LENGTH];
    uint32_t mapq[BS_GENOMES_COUNT], mapq2[BS_GENOMES_COUNT];
    long long int tlen, tlen2;
    bzero(chosen, sizeof(chosen));

    while (true) {
        for (int i = 0; i < BS_GENOMES_COUNT; i++) {
            readPenalties[i] = 0;
            if (!ReadLine(line[i], MAX_SAM_LINE_LENGTH, samFiles[i])) return;
            sscanf(line[i], "%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", qname, &flag[i], rname[i],
                   &pos[i], &mapq[i], cigar[i], rnext, pnext, &tlen, seq_string, quality_string);
            int index = ChromIndex(rname[i]);
            if (index == -1) readPenalties[i] += UNALIGNED_PENALTY;
            else {
                int refCT = 0; // if read has been aligned to last two genome
                if (i <= 2) // if read has been aligned to first thre genomes
                    refCT = 1;
                readCigar(cigar[i], pos[i] + chrom[index].chrStart - 1, seq_string, i, index, pos[i], refCT);
            }

            if (paired) {
                readPenalties[i + BS_GENOMES_COUNT] = 0;
                if (!ReadLine(line2[i], MAX_SAM_LINE_LENGTH, samFiles[i])) {
                    cerr << "Error: unexpected end of input SAM file while reading the paired end.\n";
                    return;
                }
                sscanf(line2[i], "%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", qname2, &flagTwo[i],
                       rname2[i], &pos2[i], &mapq2[i], cigar2[i], rnext2, pnext2, &tlen2, seq_string2, quality_string2);
                int index2 = ChromIndex(rname2[i]);
                if (index2 == -1) readPenalties[i + BS_GENOMES_COUNT] += UNALIGNED_PENALTY;
                else {
                    int refCT = 0; // if read has been aligned to last two genome
                    if (i <= 2) // if read has been aligned to first thre genomes
                        refCT = 1;
                    readCigar(cigar2[i], pos2[i] + chrom[index2].chrStart - 1, seq_string2, i + BS_GENOMES_COUNT,
                              index2, pos2[i], refCT);
                }
            }
            cerr << i << "\t" << readPenalties[i] << "\t" << readPenalties[i + BS_GENOMES_COUNT] << endl;
        }

        int min, min2;
        find_min_penalties(&min, &min2);
        chosen[min]++; // shows how many times a genome has been selected
        fprintf(outputFile, "%s", line[min]);
        if (paired) {
            chosen[min2]++; // shows how many times a genome has been selected
            fprintf(outputFile, "%s", line2[min2]);
        }
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
    char *samFileName = nullptr, c;
    char *annotationFile, *referenceName;
    outputFile = stdout;
    while ((c = (char) getopt_long(argc, argv, "x:s:c:p:o:h:Pm:M:123l:", long_options, &option_index)) >= 0) {
        switch (c) {
            case 'x':
                referenceName = strdup(optarg);
                break;
            case 's':
                samFileName = strdup(optarg);
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
                mismatch_penalty = (int) strtol(optarg, &strtox_temp, 10);
                go_penatly = (int) strtol(argv[optind], &strtox_temp, 10);
                ge_penalty = (int) strtol(argv[optind + 1], &strtox_temp, 10);
                optind += 2;
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
                paired_min_distance = (int) strtol(optarg, &strtox_temp, 10);
                break;
            case 'M':
                paired_max_distance = (int) strtol(optarg, &strtox_temp, 10);
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
//                limit = (int) strtol(optarg, &strtox_temp, 10);;
            default:
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
    Process(false);
    for (int j = 0; j < BS_GENOMES_COUNT; j++)
        fprintf(stderr, "Number of reads aligned to genome %d: %d\n", j, chosen[j]);

    if (outputFile != stdout) fclose(outputFile);
}


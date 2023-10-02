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
#include <map>

#include "const.h"

using namespace std;

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

struct Chrom {
    string chrName;
    unsigned long long chrStart{};
}; // for converting refrence position to chromosome based position

vector<Chrom> chrom;
vector<string> chromName;
vector<Island> islands;
map<uint64_t, int> cpg_ids;
FILE *outputFile = nullptr, *headerFile = nullptr;
bool paired = false, em_flag = false;
int chromNum, paired_min_distance = 0, paired_max_distance = 10000,
        mismatch_penalty = 0, go_penatly = 0, ge_penalty = 0,
        index_usage_count[BS_GENOMES_COUNT];
vector<int> read_min_penalty, read_min_penalty2,
        cpg_read_count, cpg_methyl_count,
        new_cpg_read_count, new_cpg_methyl_count;
char *genome, *samNames[BS_GENOMES_COUNT], *strtox_temp;


inline char to_upper(char c) {
    if (c >= 'a')
        return (char) (c - 'a' + 'A');
    return c;
}

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
    uint64_t gs = 0;
    chromNum = 0;
    fprintf(stderr, "Reading genome...\n");
    char fLineMain[1000000];
    Chrom ch;
    int cpg_counts = 0;
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
            uint64_t i = gs > 0 ? gs : gs + 1;
            if (em_flag)
                for (; i < gs + n; ++i) {
                    genome[i] = to_upper(genome[i]);
                    if (genome[i - 1] == 'C' && genome[i] == 'G') {
                        cpg_ids[i - 1] = cpg_counts++;
                        cpg_ids[i] = cpg_counts++;
                        i++;
                    }
                }
            gs += n;
        }
    }
    fclose(f);
    return 0;
}


inline int ChromIndex(const char *chr) { //}, size_t len = 100000) { use "len" to make chr16_xxx to chr16
    const char *chr_dup = nullptr, *chr_to_compare = chr;
//    if (len < strlen(chr)) {
//        chr_dup = strndup(chr, len);
//        chr_to_compare = chr_dup;
//    }
    for (unsigned int i = 0; i < chrom.size(); i++)
        if (chrom[i].chrName == chr_to_compare)
            return (int) i;
    delete chr_dup;
    return -1;
}

void find_min_penalties(const char rname[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH],
                        const char rname2[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH], const int *flag, const int *flagTwo,
                        const uint64_t *pos, const uint64_t *pos2, const double *readPenalties, int *min1, int *min2) {
    *min1 = *min2 = 0;
    double min_penalty = readPenalties[0] + (paired ? readPenalties[0 + BS_GENOMES_COUNT] : 0);
    for (int i = 0; i < BS_GENOMES_COUNT; i++) {
        if (!paired) {
            if (min_penalty > readPenalties[i]) {
                min_penalty = readPenalties[i];
                *min1 = i;
            }
        } else {
            int chr_index = ChromIndex(rname[i]); //, strchr(rname[i], '_') - rname[i]);
            if (chr_index < 0)
                continue;
            for (int j = 0; j < BS_GENOMES_COUNT; j++) {
                double new_penalty = readPenalties[i] + readPenalties[j + BS_GENOMES_COUNT];
                int chr_index2 = ChromIndex(rname2[j]); //, strchr(rname2[j], '_') - rname2[j]);
                if (chr_index2 < 0)
                    continue;
                orientation_t orien;
                if ((flag[i] & 16) == (flagTwo[j] & 16))
                    orien = ff;
                else if (((flag[i] & 16) == 0 && pos[i] < pos2[j]) ||
                         ((flag[i] & 16) != 0 && pos[i] > pos2[j]))
                    orien = fr;
                else
                    orien = rf;
                if (chr_index != chr_index2 ||
                    llabs((signed) (pos[i] - pos2[j])) < paired_min_distance ||
                    llabs((signed) (pos[i] - pos2[j])) > paired_max_distance ||
                    (orientation != all && orientation != orien))
                    new_penalty += DISCORD_PENALTY;
                if (min_penalty > new_penalty) {
                    min_penalty = new_penalty;
                    *min1 = i;
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


double calc_base_nuc_penalty(const uint64_t ref_i, const char read_ch, const int chr, const uint64_t chrPos,
                             vector<pair<int, bool>> &read_cpg_ids, const bool refCT) {
    const char read = to_upper(read_ch);
    genome[ref_i] = to_upper(genome[ref_i]);
    genome[ref_i + 1] = to_upper(genome[ref_i + 1]);
    genome[ref_i - 1] = to_upper(genome[ref_i - 1]);
    char methyl_nuc = 'G', bis_methyl_nuc = 'A', CPG_next_nuc = 'C';
    int next_nuc = -1;
    if (refCT) {
        methyl_nuc = 'C';
        bis_methyl_nuc = 'T';
        next_nuc = 1;
        CPG_next_nuc = 'G';
    }

    if ((read != bis_methyl_nuc && read != methyl_nuc) || (genome[ref_i] != methyl_nuc)) {
        if (genome[ref_i] != read)
            return mismatch_penalty;
    } else {
        bool is_methylated = read != bis_methyl_nuc;
        double methyl_ratio = 0;
        if (genome[ref_i + next_nuc] == CPG_next_nuc) {// is CpG
            bool methyl_ratio_changed = false;
            if (!cpg_ids.empty()) {
                int cpg_id = cpg_ids[ref_i];
                read_cpg_ids.emplace_back(cpg_id, is_methylated);
                if (cpg_read_count[cpg_id] > 0) {
                    methyl_ratio = 1.0 * cpg_methyl_count[cpg_id] / cpg_read_count[cpg_id];
                    methyl_ratio_changed = true;
                }
            }
            if (!methyl_ratio_changed && !isInIsland(chr, chrPos))
                methyl_ratio = 1;
        }
        return (is_methylated ? (1 - methyl_ratio) : methyl_ratio) * mismatch_penalty;
    }
    return 0;
}

// Reads Cigar sequence and call calc_base_nuc_penalty for calculating penalties for each base
double calc_read_penalty(const char *cigar, const char *seq_string, const int chr_index, const uint64_t chr_pos,
                         vector<pair<int, bool>> &read_cpg_ids, const bool refCT) {
    int cpos = 0, value = 0;
    double result = 0;
    uint64_t ref_index = chr_pos + chrom[chr_index].chrStart - 1;
    long read_index = 0;
    while (true) {
        if (!isdigit(cigar[cpos])) {
            if (value > 0) {
                if (cigar[cpos] == 'M' || cigar[cpos] == 'm') {
                    for (int j = 0; j < value; j++) {
                        result += calc_base_nuc_penalty(ref_index, seq_string[read_index], chr_index, chr_pos,
                                                        read_cpg_ids, refCT);
                        ref_index++;
                        read_index++;
                    }
                } else if (cigar[cpos] == 'D' || cigar[cpos] == 'd') {
                    ref_index += value;
                    result += ge_penalty * (value - 1) + go_penatly;
                } else if (cigar[cpos] == 'I' || cigar[cpos] == 'i') {
                    read_index += value;
                    result += ge_penalty * (value - 1) + go_penatly;
                }
            } else if (cigar[cpos] == '*')
                return UNALIGNED_PENALTY;     // maximum penalty for not aligned reads
            if (cigar[cpos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[cpos] - '0';
        }
        cpos++;
    }
    return result;
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

double read_sam_and_calc_penalty(FILE *sam_file, char *line, int *flag, char *rname,
                                 uint64_t *pos, vector<pair<int, bool>> &read_cpg_ids, const bool refCT) {
    char cigar[MAX_READ_LENGTH * 2], qname[MAX_READ_NAME_LENGTH], rnext[MAX_CHR_NAME_LENGTH], pnext[30],
            seq_string[MAX_READ_LENGTH], quality_string[MAX_READ_LENGTH];
    uint32_t mapq;
    long long tlen;
    if (!ReadLine(line, MAX_SAM_LINE_LENGTH, sam_file))
        return -1;
    sscanf(line, "%s\t%d\t%s\t%" PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n", qname, flag, rname,
           pos, &mapq, cigar, rnext, pnext, &tlen, seq_string, quality_string);
    int chr_index = ChromIndex(rname);
    if (chr_index == -1)
        return UNALIGNED_PENALTY;
    return calc_read_penalty(cigar, seq_string, chr_index, *pos, read_cpg_ids, refCT);
}

void change_read_candidate(vector<pair<int, bool>> *read_cpg_ids, int prev_candidate, int new_candidate) {
    if (prev_candidate >= 0)
        for (auto p : read_cpg_ids[prev_candidate]) {
            new_cpg_read_count[p.first]--;
            if (p.second)
                new_cpg_methyl_count[p.first]--;
        }
    for (auto p : read_cpg_ids[new_candidate]) {
        new_cpg_read_count[p.first]++;
        if (p.second)
            new_cpg_methyl_count[p.first]++;
    }
}

uint64_t Process() { // return em total changed
    char rname[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH], rname2[BS_GENOMES_COUNT][MAX_CHR_NAME_LENGTH],
            line[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH], line2[BS_GENOMES_COUNT][MAX_SAM_LINE_LENGTH];
    int flag[BS_GENOMES_COUNT], flagTwo[BS_GENOMES_COUNT];
    uint64_t pos[BS_GENOMES_COUNT], pos2[BS_GENOMES_COUNT];
    double readPenalties[2 * BS_GENOMES_COUNT];

    FILE *samFiles[BS_GENOMES_COUNT];
    for (int i = 0; i < BS_GENOMES_COUNT; i++)
        samFiles[i] = fopen(samNames[i], "r");
    bzero(index_usage_count, sizeof(index_usage_count));
    uint64_t read_count = 0, em_total_changed = 0;
    if (cpg_read_count.empty()) {
        cpg_read_count.resize(cpg_ids.size());
        cpg_methyl_count.resize(cpg_ids.size());
        new_cpg_read_count.resize(cpg_ids.size());
        new_cpg_methyl_count.resize(cpg_ids.size());
    } else {
        cpg_read_count = new_cpg_read_count;
        cpg_methyl_count = new_cpg_methyl_count;
    }
    while (true) {
        vector<pair<int, bool>> read_cpg_ids[BS_GENOMES_COUNT], read_cpg_ids2[BS_GENOMES_COUNT];
        for (int i = 0; i < BS_GENOMES_COUNT; i++) {
            readPenalties[i] = read_sam_and_calc_penalty(samFiles[i], line[i], flag + i, rname[i], pos + i,
                                                         read_cpg_ids[i], i <= 2);
            if (readPenalties[i] < 0)
                return em_total_changed;

            if (paired) {
                readPenalties[i + BS_GENOMES_COUNT] =
                        read_sam_and_calc_penalty(samFiles[i], line2[i], flagTwo + i, rname2[i], pos2 + i,
                                                  read_cpg_ids2[i], i <= 2);
                if (readPenalties[i + BS_GENOMES_COUNT] < 0) {
                    cerr << "Error: unexpected end of input SAM file while reading the paired end.\n";
                    return em_total_changed;
                }
            }
        }

        int min, min2;
        find_min_penalties(rname, rname2, flag, flagTwo, pos, pos2, readPenalties, &min, &min2);
        if (em_flag) {
            if (read_min_penalty.size() <= read_count)
                read_min_penalty.push_back(0);
            if (read_min_penalty[read_count] != min + 1) {
                em_total_changed += 1;
                change_read_candidate(read_cpg_ids, read_min_penalty[read_count] - 1, min);
                read_min_penalty[read_count] = min + 1;
            }
            if (paired) {
                if (read_min_penalty2.size() <= read_count)
                    read_min_penalty2.push_back(0);
                if (read_min_penalty2[read_count] != min2 + 1) {
                    em_total_changed += 1;
                    change_read_candidate(read_cpg_ids2, read_min_penalty2[read_count] - 1, min2);
                    read_min_penalty2[read_count] = min2 + 1;
                }
            }
        } else {
            index_usage_count[min]++; // shows how many times a genome has been selected
            fprintf(outputFile, "%s", line[min]);
            if (paired) {
                index_usage_count[min2]++; // shows how many times a genome has been selected
                fprintf(outputFile, "%s", line2[min2]);
            }
        }
        read_count++;
    }
}


int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr,
                "Usage: align_bs -x <reference genome> -s <input sam files> -c <CpG-island map> -p <penalties> [-o <output SAM file>] [-h <SAM header file>] [-l/--limit <maximum allowed mismatches>]\n");
        fprintf(stderr, "Additional arguments for paired-end reads:\n");
        fprintf(stderr,
                "    -P to consider the input as paired-end, -m/--min for minimum and -M/--max for maximum distance, --fr/rf/ff for orientation, and -e/--em for using em to improve alignment.\n");
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
            {"em",          no_argument,       nullptr, 'e'},
            {"limit",       required_argument, nullptr, 'l'}
    };
    int option_index = 0;
    char *samFileName = nullptr, *annotationFile = nullptr, *referenceName = nullptr, c;
    outputFile = stdout;
    while ((c = (char) getopt_long(argc, argv, "x:s:c:p:o:h:Pm:M:123l:e", long_options, &option_index)) >= 0) {
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
            case 'e':
                em_flag = true;
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
    while (true) {
        unsigned long long changed_count, old_changed_count = 0;
        changed_count = Process();
        fprintf(stderr, "Number of changed reads position: %llu\n", changed_count);
        if (changed_count < 10 || (old_changed_count > 0 && changed_count > (long double) 0.9 * old_changed_count)) {
            if (em_flag)
                em_flag = false;
            else
                break;
        }
        old_changed_count = changed_count;
    }
    for (int j = 0; j < BS_GENOMES_COUNT; j++)
        fprintf(stderr, "Number of reads aligned to genome %d: %d\n", j, index_usage_count[j]);

    if (outputFile != stdout) fclose(outputFile);
}
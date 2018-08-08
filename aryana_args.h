#include <stdint.h>
#include <stdio.h>

#ifndef __aryana_args
#define __aryana_args

typedef enum {
    orien_all, orien_ff, orien_fr, orien_rf
} orientation_t;

// ignore_none: computer one mismatch for any pair of different nucleotides
// ignore_CT: do not count C->T mismatches. This is for alignment of bisulfite reads
// ignore_GA: do not count G->A mismatches. This is for alignment of the PCR product of bisulfite reads.

typedef enum {
    ignore_none, ignore_CT, ignore_GA
} ignore_mismatch_t;

typedef enum {
    illumina, pacbio, nanopore
} platform_t;

typedef struct {
    int paired;
    uint64_t min_dis, max_dis;
    int pairID;
    orientation_t orientation;
    int discordant;
    FILE **index_files;
    size_t *index_files_size;
    char *reference;
    char *read_file, *read_file2;
    int single;
    int threads;
    int potents;
    int seed_length;
    double best_factor;
    int bisulfite;
    int order;
    int debug;
    int exactmatch_num;
    int report_multi, mismatch_limit;
    int out_buffer_factor;
    int match_score, mismatch_penalty, gap_open_penalty, gap_ext_penalty;
    platform_t platform;
    ignore_mismatch_t ignore;
    int tag_size;
    double indel_ratio_between_seeds;
    FILE *stdout_file;
} aryana_args;
#endif

#include <stdint.h>
#ifndef __aryana_args
#define __aryana_args
typedef struct {
    int paired;
    uint64_t min_dis,max_dis;
    uint64_t len1,len2;
    int pairID;
    char ori[2];
    int discordant;
    char *reference;
    char *read_file, *read_file1, *read_file2;
    int single;
    int threads;
    int potents;
    int seed_length;
    int best_factor;
    int bisulfite;
    int order;
    int debug;
} aryana_args;
#endif

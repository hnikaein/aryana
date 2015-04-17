#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "aryana_args.h"
#include "utils.h"
#include "bwa2.h"
#include "main.h"
#define aryana_version "0.1"
int debug = 0;
void bwa_print_sam_PG()
{
    printf("@PG\tID:bwa\tPN:bwa\tVN:%s\n", aryana_version);
}

int get_value_string(int argc, const char **argv, char * arg, char * value) {
    int i = 1;
    for(i=1; i<argc - 1; i++)
        if(strcmp(argv[i], arg) == 0) {
            strcpy(value, argv[i + 1]);
            //fprintf(stderr, "value = %s, arg = %s\n", value, arg);
            return 1;
        }
    return 0;
}


void Usage() {
    fprintf(stderr, "Program: aryana\n");
    fprintf(stderr, "Version: %s\n", aryana_version);
    fprintf(stderr, "Usage:\n\n");
    fprintf(stderr, "Creating index from reference genome:\n");
    fprintf(stderr, "aryana index [reference genome (fasta)]\n");
    fprintf(stderr, "aryana fa2bin [reference genome (fasta)]\n\n");
    fprintf(stderr, "Alignment of single end reads:\n");
    fprintf(stderr, "aryana [-x,--index <reference genome index>] [-i,--input <reads file (fastq format)>] [-o,--output <alignment file (SAM format)>]\n");
    fprintf(stderr, "Optional arguments: \n");
    fprintf(stderr, "[-t,--threads <threads number>] [-c,--candidates <num of alignment candidates, default=10>] [-s,--seed <fixed length of seed sequence, default=dynamic>]\n");
    fprintf(stderr, "[-D <debug info level, default=0>] [-O (keep the order of reads in output as the input)]\n\n");
    fprintf(stderr, "Alignment of paired end reads:\n");
    fprintf(stderr, "aryana [-x,--index <reference genome index>] [-1,--first <reads file 1 (fastq format)>] [-2,--second <reads file 2 (fastq format)>]\n");
    fprintf(stderr, "Optional arguments:\n");
    fprintf(stderr, "[--fr, --ff, --rf (orientation of paired ends)] [-m,--min <min distance between pair reads>] [-M,--max  <max distance between pair reads>]\n");
    fprintf(stderr, "[-d,--no-discordant (do not print discordants reads)]\n\n");
    fprintf(stderr, "Alignment of bisulfite-sequencing reads:\n");
    fprintf(stderr, "[-b,--bisulfite <bisulfite reference genome index>]\n\n");
    fprintf(stderr, "[-e, <number of selected exact matches>]\n\n");
    fprintf(stderr, "See README.md for more details.\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    if (argc < 3) Usage();
    aryana_args args;
    args.discordant=1;
    args.threads=1;
    args.potents=10;
    args.seed_length = -1;
    args.best_factor = -1;
    args.bisulfite = 0;
    args.order = 0;
    char *refNames[5];
    bzero(refNames, sizeof(refNames));
    if (strcmp(argv[1], "index") == 0)  return bwa_index(argc-1, argv+1);
    if (strcmp(argv[1], "fa2bin") == 0) return fa2bin(argc-1, argv+1);
    static struct option long_options[] =
    {
        {"output", required_argument, 0, 'o'},
        {"index", required_argument, 0, 'x'},
        {"input", required_argument, 0, 'i'},
        {"first", required_argument, 0, '1'},
        {"second", required_argument, 0, '2'},
        {"fr", no_argument, 0, '3'},
        {"rf", no_argument, 0, '4'},
        {"ff", no_argument, 0, '5'},
        {"min", required_argument, 0, 'm'},
        {"max", required_argument, 0, 'M'},
        {"threads", required_argument, 0, 't'},
        {"seed", required_argument, 0, 's'},
        {"candidates", required_argument, 0, 'c'},
        {"factor", required_argument, 0, 'f'},
        {"bisulfite", required_argument, 0, 'b'},
        {"order", no_argument, 0, 'O'},
        {"debug", required_argument, 0, 'D'},
        {"no-discordant", no_argument, 0, 'd'},
        {"exact-match", required_argument, 0, 'e'},
    };
    char* output = NULL;
    char* inputFolder;
    int option_index = 0;
    int c;
    while((c = getopt_long(argc, argv, "o:x:i:1:2:345m:M:t:s:c:f:b:e:OD:d", long_options, &option_index)) >= 0) {
        switch(c) {
        case 'o':
            output = strdup(optarg);
            break;
        case 'x':
            args.reference = strdup(optarg);
            break;
        case 'i':
            args.read_file = strdup(optarg);
            args.single = 1;
            args.paired = 0;
            break;
        case '1':
            args.paired = 1;
            args.read_file1 = strdup(optarg);
            break;
        case '2':
            args.paired = 1;
            args.read_file2 = strdup(optarg);
            break;
        case '3':
            strcpy(args.ori, "fr");
            break;
        case '4':
            strcpy(args.ori, "rf");
            break;
        case '5':
            strcpy(args.ori, "ff");
            break;
        case 'm':
            args.min_dis = atoi(optarg);
            break;
        case 'M':
            args.max_dis = atoi(optarg);
            break;
        case 't':
            args.threads = atoi(optarg);
            break;
        case 's':
            fprintf(stderr,"Seed: %s\n", optarg);
            args.seed_length = atoi(optarg);
            break;
        case 'c':
            args.potents = atoi(optarg);
            break;
        case 'f':
            args.best_factor = atoi(optarg);
            break;
        case 'b':
            args.bisulfite = 1;
            inputFolder = (char *) malloc(strlen(optarg)+50);
            strcpy(inputFolder, optarg);
            strcat(inputFolder, "/");
            refNames[0] = (char*) malloc(strlen(inputFolder)+50);
            strcpy(refNames[0], inputFolder);
            strcat(refNames[0], "originalGenome.fa");
            refNames[1] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[1], inputFolder);
            strcat(refNames[1], "BisulfiteGenomeIslandConsideredCT.fa");
            // refNames[2] = (char*)malloc(strlen(inputFolder)+50);
            // strcpy(refNames[2], inputFolder);
            // strcat(refNames[2], "BisulfiteGenomeContextConsideredCT.fa");
            refNames[2] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[2], inputFolder);
            strcat(refNames[2], "BisulfiteGenomeCompleteCT.fa");
            refNames[3] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[3], inputFolder);
            strcat(refNames[3], "BisulfiteGenomeIslandConsideredGA.fa");
            // refNames[5] = (char*)malloc(strlen(inputFolder)+50);
            // strcpy(refNames[5], inputFolder);
            // strcat(refNames[5], "BisulfiteGenomeContextConsideredGA.fa");
            refNames[4] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[4], inputFolder);
            strcat(refNames[4], "BisulfiteGenomeCompleteGA.fa");
            break;
        case 'e':
            args.exactmatch_num = atoi(optarg);
        case 'O':
            args.order = 1;
            break;
        case 'D':
            debug = atoi(optarg);
            break;
        case 'd':
            args.discordant=0;
            break;
        default:
            fprintf(stderr, "Invalid argument: %c\n", c);
            exit(1);
        }
    }
    if (args.threads < 1) args.threads = 1;

    if(args.bisulfite) {
        if(!output) {
            fprintf(stderr, "The ouptut name should be specified\n");
            return -1;
        }
        char *output_temp = malloc(strlen(output)+5);
        int i;
        for(i=0; i<5; i++) {
            sprintf(output_temp, "%s-%d", output, i);
            FILE *sam = fopen(output_temp, "w");
            stdout = sam;
            args.reference = refNames[i];
            bwa_aln_core2(&args);
            fclose(sam);
        }
        free(output_temp);
    }
    else {
        FILE *sam;
        if(output) {
            sam = fopen(output, "w");
            stdout = sam;
        }
        bwa_aln_core2(&args);
        if(output)
            fclose(sam);
    }

    //fprintf(stderr, "ori = %s\n", args.ori);
    //bwa_aln_single(args.reference, args.fq);
    /*	pair_opt options;
    	options.paired=0;
    	options.min_dis=0;
    	options.max_dis=0;
    	if (options.paired==0)
    		bwa_aln_core2(argv[optind], argv[optind+1],NULL,NULL, opt, &options);
    	else
    		bwa_aln_core2(argv[optind], NULL, argv[optind+1], argv[optind+2], opt, &options);*/
    return 0;
}

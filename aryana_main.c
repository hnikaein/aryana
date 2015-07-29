#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "aryana_args.h"
#include "utils.h"
#include "bwa2.h"
#include "main.h"
#define aryana_version "0.1"
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
    fprintf(stderr, "    aryana index [reference genome (fasta)]\n");
    fprintf(stderr, "    aryana fa2bin [reference genome (fasta)]\n\n");
    fprintf(stderr, "Aligning single or paired-end reads:\n");
    fprintf(stderr, "    aryana -x <reference genome index> {-i <single-end reads (fastq format)> | -1 <paired-end reads 1> -2 <paired-end reads 2>} [options]*\n\n");
    fprintf(stderr, "Optional arguments: \n");
    fprintf(stderr, "    -o/--output                   output SAM file for the alignment results, default=standard output\n");
    fprintf(stderr, "    -t/--threads <int>            the number of parallel threads, default=1\n");
    fprintf(stderr, "    -s/--seed <int>               fixed length of the seed sequence, default=dynamically selected based on the length of each read\n");
    fprintf(stderr, "    -l/--limit <int>              maximum number of mismatches allowed for each alignment, default=unlimited\n");
    fprintf(stderr, "    -O/--order                    print the reads in output with the same order in the input.\n");
    fprintf(stderr, "    -B/--buffer <int>             size of output buffer. Increase it only when using -O/--order and program issues errors. default=100000\n");
    fprintf(stderr, "    --mp <int>                    maximum mismatch penalty, default=5\n");
    fprintf(stderr, "    --go <int>                    gap open penalty, default=5\n");
    fprintf(stderr, "    --ge <int>                    gap extension penalty, default=3\n");
    fprintf(stderr, "    -r/--report-multi             report multiple alignment positions for each read\n");
    fprintf(stderr, "    -e/--exact-match              number of exact matches of each seed to check, increasing it might increase accuracy, default=50\n");
    fprintf(stderr, "    -c/--candidates <int>         number of alignment position candidates to check, increasing it might increase accuracy, default=10\n");
    fprintf(stderr, "    -D/--debug <int>              the level of printing debug info, default=0 (no debug info)\n\n");
    fprintf(stderr, "Optional arguments for paired-end alignment:\n");
    fprintf(stderr, "     --fr/--ff/--rf               relative orientation of paired ends, default: no restriction on orientation. f=forward, r=reverse.\n");
    fprintf(stderr, "                                  only one of orientation arguments might be used.\n");
    fprintf(stderr, "     -m/--min <int>               minimum distance between paired ends, default=0\n");
    fprintf(stderr, "     -M/--max <int>               maximum distance between paired ends, default=10000\n");
    fprintf(stderr, "     -d/--no-discordant           do not print discordants reads, default=if a paired alignment is not found, the best hit for each read is reported.\n\n");
    fprintf(stderr, "Alignment of bisulfite-sequencing (DNA Methylation assays) reads:\n");
    fprintf(stderr, "     aryana -b <bisulfite reference genome index>  {-i <single-end reads (fastq format)> | -1 <paired-end reads 1> -2 <paired-end reads 2>} [options]*\n\n");
    fprintf(stderr, "     Bisulfite reference genome index is created by prepare_genomes.sh, please refere to the manual for more details.\n");
    fprintf(stderr, "     All optional arguments for normal reads can be used for bisulfite-sequencing reads too.\n\n");
    fprintf(stderr, "Additional optional arguments for bisulfite-sequencing reads:\n");
    fprintf(stderr, "     --ct                        ignore C->T mismatches. While using -b this argument is automatically set based on the type of converted genome\n");
    fprintf(stderr, "     --ga                        ignore G->A mismatches. This argument is also automatically set with using -b. Either --ct or --ga can be used.\n");
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
    args.debug = 0;
    args.seed_length = -1;
    args.best_factor = -1;
    args.bisulfite = 0;
    args.order = 0;
    args.exactmatch_num = 50;
    args.report_multi = 0;
    args.mismatch_limit = -1;
    args.mismatch_penalty = 5;
    args.gap_open_penalty = 5;
    args.gap_ext_penalty = 3;
    args.out_buffer_factor = 100000;
    args.ignore = ignore_none;
    args.orientation = orien_all;
    args.min_dis = 0;
    args.max_dis = 10000;
    char *refNames[5];	// Number of bisulfite-seq reference genomes
    bzero(refNames, sizeof(refNames));
    ignore_mismatch_t ignore[5]; // We should define for each bis-Seq reference genome which type of mismatch is ignored
    if (strcmp(argv[1], "index") == 0)  return bwa_index(argc-1, argv+1);
    if (strcmp(argv[1], "fa2bin") == 0) return fa2bin(argc-1, argv+1);
    static struct option long_options[] =
    {
        {"output", required_argument, 0, 'o'},
        {"index", required_argument, 0, 'x'},
        {"input", required_argument, 0, 'i'},
        {"first", required_argument, 0, '1'},
        {"second", required_argument, 0, '2'},
        {"fr", no_argument, 0, 1},
        {"rf", no_argument, 0, 2},
        {"ff", no_argument, 0, 3},
        {"min", required_argument, 0, 'm'},
        {"max", required_argument, 0, 'M'},
        {"threads", required_argument, 0, 't'},
        {"seed", required_argument, 0, 's'},
        {"candidates", required_argument, 0, 'c'},
        {"factor", required_argument, 0, 'f'},
        {"bisulfite", required_argument, 0, 'b'},
        {"order", no_argument, 0, 'O'},
        {"buffer", required_argument, 0, 'B'},
        {"debug", required_argument, 0, 'D'},
        {"no-discordant", no_argument, 0, 'd'},
        {"exact-match", required_argument, 0, 'e'},
        {"report-multi", no_argument, 0, 'r'},
        {"limit", required_argument, 0, 'l'},
        {"ct", no_argument, 0, 4},
        {"ga", no_argument, 0, 5},
        {"mp", required_argument, 0, 6},
        {"go", required_argument, 0, 7},
        {"ge", required_argument, 0, 8},
    };
    char* output = NULL;
    char* inputFolder;
    int option_index = 0;
    int c;
    while((c = getopt_long(argc, argv, "o:x:i:1:2:345m:M:t:s:c:f:b:e:OB:D:drl:\x01\x02\x03\x04\x05\x06:\x07:\x08:", long_options, &option_index)) >= 0) {
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
            args.read_file = strdup(optarg);
            break;
        case '2':
            args.paired = 1;
            args.read_file2 = strdup(optarg);
            break;
        case 1:
            args.orientation = orien_fr;
            break;
        case 2:
            args.orientation = orien_rf;
            break;
        case 3:
            args.orientation = orien_ff;
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
            ignore[0] = ignore_none;
            refNames[1] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[1], inputFolder);
            strcat(refNames[1], "BisulfiteGenomeIslandConsideredCT.fa");
            ignore[1] = ignore_CT;
            // refNames[2] = (char*)malloc(strlen(inputFolder)+50);
            // strcpy(refNames[2], inputFolder);
            // strcat(refNames[2], "BisulfiteGenomeContextConsideredCT.fa");
            refNames[2] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[2], inputFolder);
            strcat(refNames[2], "BisulfiteGenomeCompleteCT.fa");
            ignore[2] = ignore_CT;
            refNames[3] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[3], inputFolder);
            strcat(refNames[3], "BisulfiteGenomeIslandConsideredGA.fa");
            ignore[3] = ignore_GA;
            // refNames[5] = (char*)malloc(strlen(inputFolder)+50);
            // strcpy(refNames[5], inputFolder);
            // strcat(refNames[5], "BisulfiteGenomeContextConsideredGA.fa");
            refNames[4] = (char*)malloc(strlen(inputFolder)+50);
            strcpy(refNames[4], inputFolder);
            strcat(refNames[4], "BisulfiteGenomeCompleteGA.fa");
            ignore[4] = ignore_GA;
            break;
        case 'e':
            args.exactmatch_num = atoi(optarg);
            break;
        case 'r':
            args.report_multi = 1;
            break;
        case 'O':
            args.order = 1;
            break;
        case 'B':
            args.out_buffer_factor = atoi(optarg);
            break;
        case 'D':
            args.debug = atoi(optarg);
            break;
        case 'd':
            args.discordant=0;
            break;
        case 'l':
            args.mismatch_limit = atoi(optarg);
            break;
        case 4:
            args.ignore = ignore_CT;
            break;
        case 5:
            args.ignore = ignore_GA;
            break;
        case 6:
            args.mismatch_penalty = atoi(optarg);
            break;
        case 7:
            args.gap_open_penalty = atoi(optarg);
            break;
        case 8:
            args.gap_ext_penalty = atoi(optarg);
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
            args.ignore = ignore[i];
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

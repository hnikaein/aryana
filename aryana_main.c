#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "aryana_args.h"
#include "utils.h"
#include "bwa2.h"
#include "main.h"
#define aryana_version "1.0.0"

void bwa_print_sam_PG()
{
	        printf("@PG\tID:bwa\tPN:bwa\tVN:%s\n", aryana_version);
}

int get_value_string(int argc, const char **argv, char * arg, char * value){
	int i = 1;
	for(i=1;i<argc - 1;i++)
		if(strcmp(argv[i], arg) == 0){
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
	fprintf(stderr, "aryana [-x Index] [-U input (fastq)]\n");
	fprintf(stderr, "Optional arguments: [-p INT (threads number)] [-P INT (potential candidates)] [--seed INT (seed length)]\n\n");
	fprintf(stderr, "Alignment of paired end reads:\n");
	fprintf(stderr, "aryana [-x Index] [-1 input1 (fastq)] [-2 input2 (fastq)]\n");
	fprintf(stderr, "Optional arguments: --{fr, ff, rf (orientation of paired ends)} [-I INT (min distance)] [-X INT (max distance)]\n");
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
	//fprintf(stderr,"main:\n");
	/*char *reference = (char*)calloc(1000, 1);
	if(!get_value_string(argc, argv, "-x", reference)){
		fprintf(stderr, "Reference genome not specified.\n");
		return -2;
	}
	char *fastq = (char*)calloc(1000, 1);
	if(!get_value_string(argc, argv, "-U", fastq)){
		fprintf(stderr, "Read sequencse file not specified.\n");
		return -2;
	}*/
	//fprintf(stderr, "fastq = %s\n", fastq);
	static struct option long_options[] =
	{
		{"qseq", no_argument, 0, 'Q'},//Reads are QSEQ files.
		{"skip", no_argument, 0, 's'},
		{"qupto", no_argument, 0, 'u'},
		{"trim5", no_argument, 0, '5'},
		{"trim3", no_argument, 0, '3'},
		{"phred33", no_argument, 0, 'M'},
		{"phred64", no_argument, 0, '6'},
		{"very-fast", no_argument, 0, 'V'},
		{"fast", no_argument, 0, 'F'},
		{"sensitive", no_argument, 0, 'Y'},
		{"very-sensitive", no_argument, 0, 'Z'},
		{"dpad", no_argument, 0, 'D'},
		{"gbar", no_argument, 0, 'W'},
		{"ignore-quals", no_argument, 0, 'T'},
		{"nofw", no_argument, 0, 'P'},
		{"norc", no_argument, 0, 'R'},
		{"minins", required_argument, 0, 'I'},
		{"maxins", required_argument, 0, 'X'},
		{"fr", no_argument, 0, '7'},
		{"rf", no_argument, 0, '8'},
		{"ff", no_argument, 0, '9'},
		{"no-mixed", no_argument, 0, 'G'},
		{"no-discordant", no_argument, 0, 0},
		{"dovetail", no_argument, 0, 1},
		{"no-contain", no_argument, 0, 2},
		{"no-overlap", no_argument, 0, 3},
		{"time", no_argument, 0, 't'},
		{"quiet", no_argument, 0, 4},
		{"threads", required_argument, 0, 'p'},
		{"reorder", no_argument, 0, 5},
		{"mm", no_argument, 0, 6},
		{"version", no_argument, 0, 8},
		{"help", no_argument, 0, 9},
		{"seed", required_argument, 0, 10},
		{"factor", required_argument, 0, 'F'},
		{"bisulfite", no_argument, 0, 'b'},
		{"bisulfit-refs", required_argument, 0, 'B'},
		{"output", required_argument, 0, 'o'},
		{"order", no_argument, 0, 'O'},
	//	{0, 0, 0, 0}
	};
	char* output = NULL;
	char* inputFolder;
	int option_index = 0;
	int c;
	while((c = getopt_long(argc, argv, "x:1:2:U:S:qfrcs:u:5:3:N:L:k:I:X:tp:hP:R:bB:o:O", long_options, &option_index)) >= 0){
		switch(c){
			case 0:
				args.discordant=0;
				break;
			case 1:
				break;
			case 'o':
				output = strdup(optarg);
				break;
			case 'x':
				args.reference = strdup(optarg);
				break;
			case 'U':
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
			case '7':
				strcpy(args.ori, "fr");
				break;
			case '8':
				strcpy(args.ori, "rf");
				break;
			case '9':
				strcpy(args.ori, "ff");
				break;
			case 'I':
				args.min_dis = atoi(optarg);
				break;
			case 'X':
				args.max_dis = atoi(optarg);
				break;
			case 'p':
				args.threads = atoi(optarg);
				break;
			case 'P':
				args.potents = atoi(optarg);
				break;
			case 10:
				fprintf(stderr,"Seed: %s\n", optarg);
				args.seed_length = atoi(optarg);
				break;
			case 'F':
				args.best_factor = atoi(optarg);
				break;
			case 'b':
				args.bisulfite = 1;
				break;
			case 'B':
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
			case 'O': 
				args.order = 1;
				break;
		}
	}
	if (args.threads < 1) args.threads = 1;

	if(args.bisulfite){
		if(!output){
			fprintf(stderr, "The ouptut name should be specified\n");
			return -1;
		}
		char *output_temp = malloc(strlen(output)+5);
		int i;
		for(i=0; i<5; i++){
			sprintf(output_temp, "%s-%d", output, i);
			FILE *sam = fopen(output_temp, "w");
			stdout = sam;
			args.reference = refNames[i];
			bwa_aln_core2(&args);
			fclose(sam);
		}
		free(output_temp);
	}
	else{
		FILE *sam;
		if(output){
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

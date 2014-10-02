#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

#include "utils.h"
#include "bwt.h"
//#include "align_bisulfite.h"

#include <stdio.h>
#include <stdlib.h>

#define maxGenomeSize 1e5
#define maxChromosomeNum 1000
#define numberOfGenomes 7
//const long maxGenomeSize = 4e9;
//const int maxChromosomeNum = 1000;

const int interval_side = 3000;
long chromPos[maxChromosomeNum];
unsigned long gs;
char chromName[maxChromosomeNum][100];
char * genome;
int chromNum;
char * genomeFile, *annotationFile, *outputFile;
uint64_t islandStarts[(long) maxGenomeSize];
uint64_t islandEnds[(long) maxGenomeSize];
long islandsNum = 0;
int highPenalty = 0, medPenalty = 0, lowPenalty = 0;
long penalties;
long readNum;
long readPenalties[numberOfGenomes];

char *referenceName, *annotationFile;
char *samNames[numberOfGenomes];
FILE *samFiles[numberOfGenomes];


struct Chrom
{
    char * chrName;
    long chrStart;
};

struct Chrom chrom[maxChromosomeNum] ;
char * reference;
uint32_t reference_size;
uint32_t reference_reminder;

struct ChrIslands
{
    int chrNum;
    int islandsNum;
    char * chrName;
    uint64_t islandStarts[(long) maxGenomeSize];
    uint64_t islandEnds[(long) maxGenomeSize];
    
};

struct ChrIslands chrIslands[maxChromosomeNum] ;

int count =0; ////////////////////////////////
int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr, "Need more inputs\n");
        return -1;
    }
    
    static struct option long_options[] = {
        { "samfiles", required_argument, 0, 's' },
        { "reference", required_argument, 0, 'x' },
        { "CpG-islands", required_argument, 0, 'c' },
        { "penalties", required_argument, 0, 'p' }
        //  {0, 0, 0, 0}
    };
    int option_index = 0;
	char* samFileName;
    int c, i;
    while ((c = getopt_long(argc, argv, "x:s:c:p:", long_options, &option_index))
           >= 0) {
        //printf("heu0");
        switch (c) {
            case 0:
                break;
            case 1:
                break;
            case 'x':
                referenceName = (char *) malloc(strlen(optarg));
                strcpy(referenceName, optarg);
                break;
            case 's':
                samFileName = (char *) malloc(strlen(optarg));
                strcpy(samFileName, optarg);
				for(i=0; i<numberOfGenomes; i++){
	                samNames[i] = (char *) malloc(strlen(samFileName) + 5);
	                sprintf(samNames[i], "%s-%d", samFileName, i);
				}
                break;
            case 'c':
                annotationFile = (char *) malloc(strlen(optarg));
                strcpy(annotationFile, optarg);
                break;
            case 'p':
                lowPenalty = atoi(optarg);
                medPenalty = atoi(argv[optind]);
                highPenalty = atoi(argv[optind + 1]);
                optind = optind + 2;
                break;
        }
    }
    // printf("heu1");
    ReadCpGIslands(annotationFile);
    int p=0;
    // for(p;p<100;p++)
    // fprintf(stderr,"start:  %lld    end : %lld number:  %d \n",chrIslands[4].islandStarts[p],chrIslands[4].islandEnds[p],chrIslands[4].chrNum);
    ReadGenome(referenceName);
    // printf("heu3");
    
    char * command = "sort -k1 ";
    char cmd_pointer[strlen(samNames[0]) + strlen(command) + 30];
    strcpy(cmd_pointer, command);
    strcat(cmd_pointer, samNames[0]);
    strcat(cmd_pointer, " > samFile1.sam");
    system(cmd_pointer);

    char cmd_pointer2[strlen(samNames[1]) + strlen(command) + 1];
    strcpy(cmd_pointer2, command);
    strcat(cmd_pointer2, samNames[1]);
    strcat(cmd_pointer2, " > samFile2.sam");
    system(cmd_pointer2);

    char cmd_pointer3[strlen(samNames[2]) + strlen(command) + 1];
    strcpy(cmd_pointer3, command);
    strcat(cmd_pointer3, samNames[2]);
    strcat(cmd_pointer3, " > samFile3.sam");
    system(cmd_pointer3);
    
    char cmd_pointer4[strlen(samNames[3]) + strlen(command) + 1];
    strcpy(cmd_pointer4, command);
    strcat(cmd_pointer4, samNames[3]);
    strcat(cmd_pointer4, " > samFile4.sam");
    system(cmd_pointer4);
    
    char cmd_pointer5[strlen(samNames[4]) + strlen(command) + 1];
    strcpy(cmd_pointer5, command);
    strcat(cmd_pointer5, samNames[4]);
    strcat(cmd_pointer5, " > samFile5.sam");
    system(cmd_pointer5);
    
    char cmd_pointer6[strlen(samNames[5]) + strlen(command) + 1];
    strcpy(cmd_pointer6, command);
    strcat(cmd_pointer6, samNames[5]);
    strcat(cmd_pointer6, " > samFile6.sam");
    system(cmd_pointer6);
    
    char cmd_pointer7[strlen(samNames[6]) + strlen(command) + 1];
    strcpy(cmd_pointer7, command);
    strcat(cmd_pointer7, samNames[6]);
    strcat(cmd_pointer7, " > samFile7.sam");
    system(cmd_pointer7);
    

    samFiles[0] = fopen("samFile1.sam", "r");
    samFiles[1] = fopen("samFile2.sam", "r");
    samFiles[2] = fopen("samFile3.sam", "r");
    samFiles[3] = fopen("samFile4.sam", "r");
	samFiles[4] = fopen("samFile5.sam", "r");
	samFiles[5] = fopen("samFile6.sam", "r");
	samFiles[6] = fopen("samFile7.sam", "r");
    
//    samFiles[0] = fopen(samNames[0], "r");
//    samFiles[1] = fopen(samNames[1], "r");
//    samFiles[2] = fopen(samNames[2], "r");
//    samFiles[3] = fopen(samNames[3], "r");
//    samFiles[4] = fopen(samNames[4], "r");
//    samFiles[5] = fopen(samNames[5], "r");
//    samFiles[6] = fopen(samNames[6], "r");
    
    char line[1000];
    int header = 1;
    char *qname, *rnext, *pnext, *seq_string, *quality_string;
    char *rname[numberOfGenomes], *cigar[numberOfGenomes];
    int flag[numberOfGenomes];
    uint64_t pos[numberOfGenomes];
    uint32_t mapq[numberOfGenomes];
    long long int tlen;
    for (i = 0; i < numberOfGenomes; i++) {
        rname[i] = malloc(100 * sizeof(char));
        cigar[i] = malloc(200 * sizeof(char));
    }
    qname = malloc(100 * sizeof(char));
    rnext = malloc(100 * sizeof(char));
    pnext = malloc(100 * sizeof(char));
    seq_string = malloc(1000 * sizeof(char));
    quality_string = malloc(500 * sizeof(char));
    int chosen[numberOfGenomes];
    int j=0;
    for (j; j < numberOfGenomes; j++)
		chosen[j] = 0;

	int stop = 0;

	while (1 && !stop) {
		for (i = 0; i < numberOfGenomes && !stop; i++) {
			if (fgets(line, 1000, samFiles[i]) == NULL) {
				stop = 1;
				break;
			}
			while (line[0] == '@'){
				if(fgets(line, 1000, samFiles[i]) == NULL){
					stop = 1;
					break;
				}
			}
			if(stop)
				break;
			sscanf(line,"%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag[i], rname[i], &pos[i],&mapq[i], cigar[i],rnext,pnext, &tlen,seq_string,quality_string);
			int index = ChromIndex(rname[i]);
			//fprintf(stderr, "INDEX %d\n", index);
			if(index == -1){
				readPenalties[i] += LONG_MAX;
				continue;
			}
			//fprintf(stderr, "AAA %s\n", qname);
            //printf("cigar : %s \n",cigar[i]);
            int flag2 = 0;
//             if((i <=3 && flag[i] == 0) || (flag[i] == 16 && i > 3))
//                 flag2 = 1;
//             if((i > 3 && flag[i] == 0) || (flag[i] == 16 && i <= 3))
//                 flag2 = 0;
            if(i <=3)
                flag2 = 1;
            if(i > 3)
                flag2 = 0;

            if(flag[i] == 16)
                reverseRead(seq_string);
			readCigar(cigar[i], pos[i]+chrom[index].chrStart-1, seq_string, i ,rname[i],pos[i],flag[i],flag2);
		}
		if(stop)
			break;
		int min = min_penalty();
		chosen[min]++;
		fprintf(stdout, "%s\t%d\t%s\t%"PRIu64"\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\t%d\n",qname, flag[min], rname[min], pos[min],mapq[min], cigar[min],rnext,pnext, tlen,seq_string,quality_string, min);
        int j=0;
        for (j; j < numberOfGenomes; j++)
            readPenalties[j]=0;
    }
    
    for (j=0; j < numberOfGenomes; j++)
		fprintf(stderr, "CHOSEN %d: %d\n", j, chosen[j]);

}

int ReadGenome(char * genomeFile) {
    fprintf(stderr, "Allocating memory...\n");
    struct stat file_info;
    if (stat(genomeFile, &file_info) == -1) {
        fprintf(stderr,
                "Could not get the information of file %s\nplease make sure the file exists\n",
                genomeFile);
        return -1;
    }
    FILE *fp;
    fp = fopen(genomeFile, "r");
    off_t file_size_bytes = file_info.st_size;
    reference_size = ceil(((double) file_size_bytes) / (double) (sizeof(char)));
    reference = (char *) malloc(reference_size * sizeof(char));
    memset(reference, 0, reference_size * sizeof(char));
	gs = 0;
	chromNum = 0;
    //fprintf(stderr, "Reading genome...\n");
    char fLine[10000];
    while (! feof(fp)) {
        // fprintf(stderr, "Reading genome1...\n");
		int n = fscanf(fp, "%s\n", fLine);
		if (n == EOF) break;
		n = strlen(fLine);
        // fprintf(stderr, "Reading genome2 :n %d ...\n",n);
		if (fLine[0] == '>') {
			//chromPos[chromNum++] = gs;
            chrom[chromNum++].chrStart = gs;
            char * temp = fLine;
            temp++;
			int lenght = strlen(temp);
			chrom[chromNum-1].chrName = (char *) malloc(lenght * sizeof(char));
            strcpy(chrom[chromNum-1].chrName, temp);
            //fprintf(stderr, "Reading genome3...\n");
			if (chromNum >= 1)
			    fprintf(stderr, "chrName: %s, chrStart: %ld\n",chrom[chromNum-1].chrName, chrom[chromNum-1].chrStart);
			fprintf(stderr, "%s\n",fLine);
			//strcpy(chromName[chromNum - 1], fLine);
		} else {
            //            ToUpper(fLine);
			memcpy(reference+gs, fLine, n);
			gs += n;
		}
	}
    chromPos[chromNum] = gs;
    fprintf(stderr, "Lenght: %ld\n",chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(fp);
}

void reverseRead(char *s){////
    char  c;
    char* s0 = s - 1;
    char* s1 = s;
    while (*s1) ++s1;
    while (s1-- > ++s0)
    {
        c   = *s0;
        *s0 = *s1;
        *s1 =  c;
    }
    
    return;
}

int min_penalty(){
    int i = 0,j;
    if(count++<50)
        for(j=0;j<4;j++)
            fprintf(stderr, "Penalties %d : %lld \n",j,readPenalties[j]);
    long min = readPenalties[0];
    if(min > readPenalties[1]){
        min = readPenalties[1];
        i = 1;
    }
    if(min > readPenalties[2]){
        min = readPenalties[2];
        i = 2;
    }
    if(min > readPenalties[3]){
        min = readPenalties[3];
        i = 3;
    }
    if(min > readPenalties[4]){
        min = readPenalties[4];
        i = 4;
    }
    if(min > readPenalties[5]){
        min = readPenalties[5];
        i = 5;
    }
    if(min > readPenalties[6]){
        min = readPenalties[6];
        i = 6;
    }
    return i;
}

char getNuc(uint64_t place,int flag) {
	// char atomic[4] = { 'A', 'C', 'G', 'T' };
// 	int rev = 0;
// 	fprintf(stderr, "AAAA %" PRIu64 " %" PRIu64 "\n", place, seq_len);
// 	if(place > (seq_len / 2))
// 	{
// 		place = (seq_len / 2) - (place - (seq_len / 2))-1;
// 		rev=1;
// 	}
// 	uint64_t block = place / (sizeof(bwtint_t) * 4);
// 	fprintf(stderr, "BBBB %" PRIu64 "\n", block);
// 	int offset = place % (sizeof(bwtint_t) * 4);
// 	fprintf(stderr, "CCCC %d\n", offset);
// 	uint64_t mask = 3;
// 	mask = mask & (reference[block] >> (2 * offset));
// 	fprintf(stderr, "DDDD %" PRIu64 "\n", mask);
// 	if (rev == 1)
// 		mask = 3 - mask;
// 	return atomic[mask];
    char temp;
    if (flag == 16){
        switch (reference[place]) {
            case 'T':
                return 'A';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'A':
                return 'G';
                
            default:
                break;
        }
        
    }else
        
	return reference[place];
}

inline void ToLower(char * s) {
    int i = 0;
    for (i = 0; i < strlen(s); i++)
        if (s[i] >= 'A' && s[i] <= 'Z')
            s[i] += 'a' - 'A';
}

int ChromIndex(char * chr) {
	int i = 0;
	for(i; i < maxChromosomeNum; i++){
		if(chrom[i].chrName == NULL)
			break;
		if(strcmp(chrom[i].chrName, chr) == 0)
			return i;
	}
	return -1;

}

void ReadCpGIslands(char * annotationFile) {
    // cerr << "Processing CpG island locations from file: " <<  annotationFile << endl;
    // ifstream f(annotationFile);
    //    if (! f.is_open()) {
    //       // cerr << "Error: CpG island locations file not found or could not be opened" << endl;
    //        exit(1);
    //    }
    fprintf(stderr, "salam");
    FILE* file = fopen(annotationFile, "r"); /* should check the result */
    fprintf(stderr, "salam2");
    char fLine[10000], chrom[10], strand[10];
    uint64_t wStart, wEnd, chr;
    int bin;
    
    long index = 0;
    fgets(fLine, sizeof(fLine), file);
    fgets(fLine, sizeof(fLine), file);
    wStart = 0;
    sscanf(fLine, "%s %" PRIu64 " %" PRIu64 " ", chrom, &wStart, &wEnd);
    fprintf(stderr, "%s \n",fLine);
    long chrIndex =0;
    int stop = 0;
    while (!stop) {
        
        chrIslands[chrIndex].chrName = (char *) malloc(strlen(chrom) * sizeof(char));
        strcpy(chrIslands[chrIndex].chrName, chrom);
        
        char *copy = (char *)malloc(strlen(chrom));
        memcpy(copy, chrom,strlen(chrom));

        index = 0;
        while(!strcmp(copy, chrom)){
            
            if (!wStart)
                continue;
            chrIslands[chrIndex].islandStarts[index] = wStart;
            chrIslands[chrIndex].islandEnds[index] = wEnd;
            //printf("starts: %" PRIu64 "     ends : %" PRIu64 "", chrIslands[chrIndex].islandStarts[index],chrIslands[chrIndex].islandEnds[index]);
            chrIslands[chrIndex].chrNum = ChromIndex(chrom);
            chrIslands[chrIndex].islandsNum++;
            index++;
            if(NULL==fgets(fLine, sizeof(fLine), file)){
                stop =1;
                break;
            }
            sscanf(fLine, "%s %" PRIu64 " %" PRIu64 " ", chrom, &wStart, &wEnd);
            
        }
        chrIndex++;
        index=0;
        
    }
    
//    int k=0,m=0;
//    for(k ; k< 100;k++){
//        fprintf(stderr, "number of %s: %d\n",chrIslands[k].chrName,chrIslands[k].islandsNum);
////        for(m ; m< chrIslands[k].islandsNum;m++)
////
////       fprintf(stderr,"starts: %" PRIu64 "     ends : %" PRIu64 "", chrIslands[k].islandStarts[m],chrIslands[k].islandEnds[m]);
//
//    }
}

void setPenalties(int p1, int p2, int p3) {
    highPenalty = p1;
    medPenalty = p2;
    lowPenalty = p3;
    return;
}

int isInIsland(uint64_t ref_i , char *chr) {
    
    //printf("island refindex : %" PRIu64 "\n",ref_i);
    int i=0;
    int chr2;
    for(i; i < maxChromosomeNum; i++){
		if(chrIslands[i].chrName == NULL)
			break;
		if(strcmp(chrIslands[i].chrName, chr) == 0)
			chr2=i;
	}
    if(chr2 == -1){
        fprintf(stderr, "wrong chromosome..\n");
        return -1;
    }
    uint64_t first = 0, last = chrIslands[chr2].islandsNum - 1;
    uint64_t middle = (first + last) / 2;
    int isInIsland = 0;
    
//    while (first <= last) {
//        //printf("first :%" PRIu64 "  last :%" PRIu64 "",islandStarts[middle],islandEnds[middle]);
//        if(count <2000)
//            fprintf(stderr,"chr : %d chr: %s, middle: %d",chr2,chr,middle);
//        if (chrIslands[chr2].islandStarts[middle] <= ref_i && ref_i <= chrIslands[chr2].islandEnds[middle]) {
//            
//            isInIsland = 1;
//            return 1;
//        } else if (chrIslands[chr2].islandEnds[middle] < ref_i)
//            first = middle + 1;
//        else
//            last = middle - 1;
//        
//        middle = (first + last) / 2;
//    }
    while(first < chrIslands[chr2].islandsNum){
        if (chrIslands[chr2].islandStarts[first] <= ref_i && ref_i <= chrIslands[chr2].islandEnds[first]) {
        
                        isInIsland = 1;
                        break;
        }
        first++;

    }
    if (count < 10000) {
        fprintf(stderr, "%s ref_i:%" PRIu64 " isinisland: %d\n ",chr,ref_i,isInIsland);
    }
    return isInIsland;
}


void CalcPenalties(uint64_t ref_i, char read, long readNum,char *chr,uint64_t chrPos,int flag,int flag2) {
	//printf("   7salam");
	//printf("read : %c   refrence: %c \n ", read,getNuc(ref_i, seq_len));
	//printf("read : %c   ", read);
    read = toupper(read);
    reference[ref_i] = toupper(reference[ref_i]);
    reference[ref_i+1] = toupper(reference[ref_i+1]);
    
    
	char atomic[4] = { 'A', 'C', 'G', 'T' };
    if(count++<100)
        fprintf(stderr,"read : %c    ref : %c  refindex : %" PRIu64 "  %s \n",read,getNuc(ref_i+1, flag) ,ref_i,chr);
    if(flag2){
        if (read == 'A' || read == 'G') 
            if (getNuc(ref_i, flag) != read)
                readPenalties[readNum] += highPenalty;
        else { //read = C or T

            if (read == 'T' && getNuc(ref_i, flag) == 'C') {
                if (getNuc(ref_i+1, flag) == 'G') { // in the CpG context
                    if (isInIsland(chrPos,chr)) // in CpG and also island
                        readPenalties[readNum] += medPenalty;
                    else
                        readPenalties[readNum] += highPenalty;
                
                } else // out of CpG context

                    readPenalties[readNum] += lowPenalty;
        
            } else if (read == 'C' && getNuc(ref_i, flag) == 'C') {
                if (getNuc(ref_i+1, flag) == 'G') { // in the CpG context
                    int temp = isInIsland(chrPos,chr);
                    if (temp == 1)  // in CpG and also island
                        readPenalties[readNum] += medPenalty;
                    else
                        readPenalties[readNum] += lowPenalty;
                    
                } else // out of CpG context
                    readPenalties[readNum] += highPenalty;

            } else if (read == 'C' && getNuc(ref_i+1, flag) == 'T')
                readPenalties[readNum] += highPenalty;

        }
    }
    else{
        if (read == 'C' || read == 'T') 
            if (getNuc(ref_i, flag) != read)
                readPenalties[readNum] += highPenalty;
        else { //read = C or T
            
            if (read == 'A' && getNuc(ref_i, flag) == 'G') {
                if (getNuc(ref_i+1, flag) == 'C') { // in the CpG context
                    if (isInIsland(chrPos,chr)) { // in CpG and also island
                        readPenalties[readNum] += medPenalty;
                    } else {
                        readPenalties[readNum] += highPenalty;
                    }
                } else { // out of CpG context
                    
                    readPenalties[readNum] += lowPenalty;
                }
            } else if (read == 'G' && getNuc(ref_i, flag) == 'G') {
                if (getNuc(ref_i+1, flag) == 'C') { // in the CpG context
                    int temp = isInIsland(chrPos,chr);
                    if (temp == 1) { // in CpG and also island
                        readPenalties[readNum] += medPenalty;
                    } else {
                        readPenalties[readNum] += lowPenalty;
                    }
                } else { // out of CpG context
                    readPenalties[readNum] += highPenalty;
                }
                
            } else if (read == 'G' && getNuc(ref_i+1, flag) == 'A')
                readPenalties[readNum] += highPenalty;
        }

    }

}
void readCigar(char * cigar, uint64_t ref_i, char *seq_string, long readNum,char *chr,uint64_t chrPos,int flag, int refNum) {
    //fprintf(stderr, "salam\n");
    int pos = 0;
    int value = 0;
    uint64_t ref_index = ref_i;
    long read_index = 0;
    char alignType;
    //printf("   %s\n", seq_string);
    //printf("cigar:   %s\n", cigar);
	while (1) {
		if (!isdigit(cigar[pos])) {
			//printf("   1salam\n");
			if (value > 0) {

				//printf("value:   %d",value);
				if (cigar[pos] == 'm') {
					int j;
					for (j = 0; j < value; j++) {

						CalcPenalties(ref_index, seq_string[read_index], readNum,chr,chrPos,flag,refNum);
                        ref_index++;
                        read_index++;

                        //          if(strstr(cigar,"79m1d9m1d5m1i6m"))
                        //                  fprintf(stderr,"   penalties for cigar : %lld   ",readPenalties[readNum]);
                    }
                } else if (cigar[pos] == 'd') {
                    ref_index += value;
                    readPenalties[readNum] += highPenalty * value;
                    
                } else if (cigar[pos] == 'i'){
                    read_index += value;
                    readPenalties[readNum] += highPenalty * value;
                    
                }
                
            }
            else if(cigar[pos]=='*'){
                readPenalties[readNum] += LONG_MAX;
                break;
            }
            //          printf("*");
            if (cigar[pos] == 0)
                break;
            //            else{
            //                printf("   3salam");
            //                alignType = cigar[pos];
            //                printf("  %c",alignType);
            //            }
            value = 0;
        } else {
            value = value * 10 + cigar[pos] - '0';
        }
        pos++;
    }
    if(count <20)
        fprintf(stderr, "read : %s \n, cigar : %s \n , penalties : %ld \n",seq_string,cigar,readNum);
}

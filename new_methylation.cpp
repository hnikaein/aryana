//
//  methylation.cpp
//  aryana
//
//  Created by Maryam Rabiee on 9/20/14.
//  Copyright (c) 2014 Maryam Rabiee. All rights reserved.
//


#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <math.h>

#include <iostream>
#include <fstream>
#include "methylation.h"

using namespace std;



#define maxGenomeSize 1e5
#define maxChromosomeNum 1000
#define READ_LENGHT 100
#define READ_SAM_BY 10
#define BUFFER_SIZE 200
unsigned long gs;
int chromNum;
int debug=0;

char * reference;
uint32_t reference_size;
uint32_t reference_reminder;
FILE *samFile;

void log(char* s){
#ifdef DEBUG
    cout << s << endl;
#endif
}

struct Line{
    char line[1000];
    char chr[30];
    uint64_t start , pos;
    // char *qname, *seq_string,*rname, *cigar;
    char cigar[300];
    char cigar2[500];
    char rname[200];
    string seq_string;
    char strand;
    
};
//struct line lines[1000];
std::vector<Line> lines;
//std::vector<Line> lines_neg;
struct Chrom
{
    char * chrName;
    long chrStart;
    long chrNum;
};

struct Chrom chrom[maxChromosomeNum] ;


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

int checkGAorCT2(Line line,int flag){
    int GA =0;
    int CT = 0 , i = 0;
    int refPos , readPos = 0;
    long ref = chrom[ChromIndex(line.chr)].chrStart +line.pos;
    int temp;
    while (i < line.seq_string.size()) {
        if(flag==16)
            temp = line.seq_string.size()-i-2;
        else
            temp = i -1;
        line.seq_string[i]=toupper(line.seq_string[i]);
        
        //        if(debug)
        //            cerr<<"  "<<line.seq_string[i]<<i<<reference[ref+temp]<<" ";
        if(line.seq_string[i] == 'A'){
            if(toupper(reference[ref+temp]) == 'G' && flag ==0 )
                GA++;
            else if(toupper(reference[ref+temp]) == 'C' && flag ==16 )
                GA++;
        }
        if(line.seq_string[i] == 'T'){
            if(toupper(reference[ref+temp]) == 'C'  && flag ==0)
                CT++;
            else if(toupper(reference[ref+temp]) == 'G' && flag ==16 )
                CT++;
        }
        i++;
        
    }
    if(debug){
        cerr<<"GA  :"<<GA<<endl;
        cerr<<"CT   :"<<CT<<endl;
    }
    
    
    if(GA > CT)
        return 16;
    else if(GA < CT)//////////////
        return 4;
    else
        return -1;
    
}

int count_methyl[BUFFER_SIZE][2];
long mode[200];
void ReadMethylation(Line line,bool chr_changed){
    //cerr<<"111"<<endl;
    char methyl , unmethyl;
    long refPos = chrom[ChromIndex(line.chr)].chrStart + line.pos-1;
    
    if(line.strand == '+'){
        methyl='C';
        unmethyl='T';
    }
    else{
        methyl = 'G';
        unmethyl = 'A';
    }
    cerr<<line.rname<<endl;
    for(int i=0;i < line.seq_string.size() ; i++){
        
        if(toupper(reference[refPos+i]) == methyl){
            //int index = find_position(line.pos + i, line.chr);
            //cerr<<reference[refPos+i]<<"   "<<line.seq_string[i]<<"  ";
            
            long pos = refPos+i;
            //cerr<<"aa  "<<pos<<endl;
            int index = (pos)%BUFFER_SIZE;
            cerr << "IND:"<<index<<"  ";
            if(mode[index] == -1){
                mode[index] = pos + 1 ;
            }
            else if(mode[index] != pos +1){
                cerr<<"aa  "<<pos<<endl;
                float methylation_ratio;
                methylation_ratio = ((float)count_methyl[index][0]/(count_methyl[index][0]+count_methyl[index][1]))*100.0 ;
                if((count_methyl[index][0]+count_methyl[index][1])==0)
                    methylation_ratio = 0;
                fprintf(stdout, "%s\t%ld\t%d\t%3f\n",line.chr, mode[index], count_methyl[index][0] ,methylation_ratio);
                count_methyl[index][0] = 0;
                count_methyl[index][1] = 0;
                mode[index] = pos +1;
            }
             //cerr<<"aa  "<<pos<<endl;
            if (toupper(line.seq_string[i]) == methyl)
                count_methyl[index][0]++;
            else if(toupper(line.seq_string[i]) == unmethyl)
                count_methyl[index][1]++;
            //cerr<<"bb  "<<pos<<endl;
        }
    }
}

char chromosome[30];
void compute_methylation(){
    cerr << "compute methylaaa222"<<endl;
    int result = readSamFile(samFile);
    cerr << "compute methylaaa2226"<<endl;
    bool chr_changed = false;
    //char *chr2 = new char[50];
    char chr2[50];
    cerr << lines.size();
    strcpy(chr2, lines[lines.size()-1].chr);
    strcpy(chromosome, lines[lines.size()-1].chr);
    cerr << "compute methylaaa333"<<endl;
    for (int i=0; i< lines.size(); i++) {
        if (!strcmp(chr2, lines[i].chr)) {
            strcpy(chr2, lines[i].chr);
        }
        //cerr << "compute methylaaa"<<endl;
        ReadMethylation(lines[0],chr_changed);
        lines.erase(lines.begin());
    }
    while(result != -1){
        cerr << "compute methylaaa111"<<endl;
        result = readSamFile(samFile);
        for (int i=0; i< lines.size(); i++) {
            //cerr << "compute methylaaa"<<endl;
            ReadMethylation(lines[0],chr_changed);
            lines.erase(lines.begin());
        }
    }
}

void convertRead(Line &line){
    
    int c=0;
    int j = 0;
    for(int i=0;i< line.seq_string.size();i++){
        if (line.cigar2[i] == 'i') {
            line.seq_string.erase(j,1);
            
        }
        else if (line.cigar2[i] == 'd') {
            line.seq_string.insert(j,"M");
            j++;
        }
        else
            j++;
        
    }
    
}

void reverseRead(Line &line){////
    //    if(debug)
    //        cerr <<"before :  "<<line.seq_string<<endl;
    int i,j=0 ;
    char * copy = new char[line.seq_string.size()];
    line.seq_string.copy(copy, line.seq_string.size(),0);
    for (i = line.seq_string.size()-1; i >= 0 ; i--){
        switch (copy[i]) {
            case 'A':
                line.seq_string[j]='T';
                break;
            case 'T':
                line.seq_string[j]='A';
                break;
            case 'C':
                line.seq_string[j]='G';
                break;
            case 'G':
                line.seq_string[j]='C';
                break;
                
            default:
                break;
        }
        j++;
    }
    
}


char * convertCigar(char * cigar){
    char * cigar2 = new char[500];
    int pos = 0;
    int value = 0, j , index =0;;
    long read_index = 0;
    char alignType;
    while (1) {
		if (!isdigit(cigar[pos])) {
			if (value > 0) {
                
				if (cigar[pos] == 'm')
					for (j = 0; j < value; j++) {
                        cigar2[index] = 'm';
                        index++;
                    }
                else if (cigar[pos] == 'd')
                    for (j = 0; j < value; j++) {
                        cigar2[index] = 'd';
                        index++;
                    }
                else if (cigar[pos] == 'i')
                    for (j = 0; j < value; j++) {
                        cigar2[index] = 'i';
                        index++;
                    }
                
                
            }
            if (cigar[pos] == 0)
                break;
            value = 0;
        } else {
            value = value * 10 + cigar[pos] - '0';
        }
        pos++;
    }
    return cigar2;
}

int readSamFile(FILE * samFile){
    //cerr<<"hey"<<endl;
    char line[1000];
    int header = 1;
    char *qname, *rnext, *pnext, *seq_string, *quality_string,*rname, *cigar;
    int flag, i;
    uint64_t pos;
    uint32_t mapq;
    long long int tlen;
    
    rname = new char[100];
    cigar = new char[200];
    qname = new char[100];
    rnext = new char[100];
    pnext = new char[100];
    seq_string = new char[1000];
    quality_string = new char[500];
    //cerr<<"hey"<<endl;
    int count_to=0;
    int stop = 0;
    while (count_to< READ_SAM_BY) {
        //cerr<<"hey"<<endl;
        if (fgets(line, 1000, samFile) == NULL){
            stop = -1;
            break;
        }
        
        while (line[0] == '@'){
            if(fgets(line, 1000, samFile) == NULL){
                stop = 1;
                break;
            }
        }
        //cout<<line<<endl;
        if(stop)
            break;
        struct Line temp;
        
        sscanf(line,"%s\t%d\t%s\t%lld\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",temp.rname, &flag, temp.chr, &(temp.pos) ,&mapq, temp.cigar,rnext,pnext, &tlen,seq_string,quality_string);
        string str(seq_string);
        temp.seq_string = str;
        
        char * cigar2 = new char[500];
        
        if(flag != 4){
            cigar2 = convertCigar(temp.cigar);
            strcpy(temp.cigar2, cigar2);
            if(debug)
                cerr << "hey  "<<temp.cigar2<<endl;
        }
        else{
            continue;
        }
        if(debug)
            cerr<<line<<"line:      2"<<endl;
        convertRead(temp);
        cerr<<"flag :"<<flag<<endl;
        int result = checkGAorCT2(temp,flag);
        if(debug)
            cerr << result<<endl;
        
        if (result == -1) {
            continue;//////////////////////////////////////////
        }
        if(result != 4){
            if (flag==16){
                reverseRead(temp);/////////////////////////////////////////////////////////////////
                temp.strand = '+';
            }
            else
                temp.strand = '-';
        }
        else{
            if (flag==16){
                reverseRead(temp);
                temp.strand = '-';
            }
            else
                temp.strand = '+';
        }
        cerr << "\n"<<temp.strand <<endl<<endl;
        count_to++;
        
        lines.push_back(temp);
    }//while
    if(stop == -1)
        return -1;
    return 1;
    
}

//void readSam_Compute(){
//    
//    int result = readSamFile(samFile);
//    //    for (int l = 0; l<lines.size(); l++) {
//    //        cout << lines[l].seq_string <<"   "<<lines[l].strand<<endl;
//    //    }
//    if(debug)
//        cerr<<result <<"llll        "<<lines.size()<<endl;
//    while (result == 1) {
//        computeMethylation();
//        result = readSamFile(samFile);
//        if(result == -1)
//            result = 2;
//    }
//    if(result != 2)
//        computeMethylation();
//}


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
		int n = fscanf(fp, "%s\n", fLine);
		if (n == EOF) break;
		n = strlen(fLine);
		if (fLine[0] == '>') {
			//chromPos[chromNum++] = gs;
            chrom[chromNum++].chrStart = gs;
            char * temp = fLine;
            temp++;
			int lenght = strlen(temp);
			chrom[chromNum-1].chrName = (char *) malloc(lenght * sizeof(char));
            strcpy(chrom[chromNum-1].chrName, temp);
			if (chromNum >= 1)
			    fprintf(stderr, "chrName: %s, chrStart: %ld\n",chrom[chromNum-1].chrName, chrom[chromNum-1].chrStart);
			fprintf(stderr, "%s\n",fLine);
		} else {
            //            ToUpper(fLine);
			memcpy(reference+gs, fLine, n);
			gs += n;
		}
	}
    //chromPos[chromNum] = gs;
    // fprintf(stderr, "Lenght: %ld\n",chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(fp);
}


int main(int argc, char *argv[]) {
    char *referenceName, *annotationFile;
    
    if ( argc != 4 ){
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s filename", argv[0] );
    }
    else
    {
        ReadGenome(argv[1]);
        cerr<<"complete reading genome..."<<endl;
        //        if(debug)
        //        for(int k = 0 ; k<200;k++)
        //            cerr<<reference[k]<<" ";
        samFile = fopen( argv[2], "r" );
        if ( samFile == 0 )
            printf( "Could not open file\n" );
        //
        debug = atoi(argv[3]);
        //        readSam_Compute();
        
        memset(mode, -1, sizeof(int)*BUFFER_SIZE);
        memset(count_methyl, 0, sizeof(count_methyl[0][0]) * 2 * BUFFER_SIZE);
        
        compute_methylation();

        float methylation_ratio;
        for(int i=0;i<BUFFER_SIZE ; i++){
            if(count_methyl[i][0] != 0 || count_methyl[i][1] !=0){
                
                methylation_ratio = ((float)count_methyl[i][0]/(count_methyl[i][0]+count_methyl[i][1]))*100.0 ;
                if((count_methyl[i][0]+count_methyl[i][1])==0)
                    methylation_ratio = 0;
                fprintf(stdout, "%s\t%ld\t%d\t%3f\n", chromosome,mode[i], count_methyl[i][0] ,methylation_ratio);
                
                
            }
        }
    }
    
}//main


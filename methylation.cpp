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
#define READ_SAM_BY 5000
unsigned long gs;
char chromName[maxChromosomeNum][100];
int chromNum;
int debug=0;

char * reference;
uint32_t reference_size;
uint32_t reference_reminder;
FILE *samFile;
int firstPointer = 0 , secondPointer=0;

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
struct Chrom
{
    char * chrName;
    long chrStart;
    long chrNum;
};

struct Cytosine
{
    long pos;
    char chr[30];
    int methylated;
    int unmethylated ;
    Cytosine(long p,char * chrom) {
        methylated=0;
        unmethylated = 0;
        pos = p;
        strcpy(chr, chrom);
	}
    Cytosine(){
        
    }
    
};
std::vector<Cytosine> cytosines;
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

int checkGAorCT2(Line line){
    int GA =0;
    int CT = 0 , i = 0;
    int refPos , readPos = 0;
    //cerr<<"GA5 "<<endl;
    long ref = chrom[ChromIndex(line.chr)].chrStart +line.pos;
    //cerr << "cigar2"<<line.cigar2<<endl;

    while (i < line.seq_string.size()) {
        line.seq_string[i]=toupper(line.seq_string[i]);
        if(debug)
        cerr<<"  "<<line.seq_string[i]<<i<<reference[ref+i]<<" ";
        if(line.seq_string[i] == 'A'){
            if(toupper(reference[ref+i]) == 'G')
                GA++;
        }
        if(line.seq_string[i] == 'T'){
            if(toupper(reference[ref+i]) == 'C')
                CT++;
        }
        i++;
        
    }
    if(debug){
    cerr<<"cigar2[i]  :"<<line.cigar2[i]<<" "<<i<<endl;
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


void computeMethylation(){


    if(lines.size() == 0){
        return;
    }
    int lastchecked = lines[0].pos;

    if(cytosines.size()!=0){
        if(debug){
            cerr << "heyyyy888   "<<strcmp(cytosines[cytosines.size()-1].chr,lines[0].chr)<<endl;
            cerr << lines[0].pos <<"    "<<lines[0].seq_string <<"   "<<lines[0].chr<<"  "<<cytosines[cytosines.size()-1].chr<<endl;
        }
        lastchecked= cytosines[cytosines.size()-1].pos;
        if (lines[0].pos > lastchecked || strcmp(cytosines[cytosines.size()-1].chr,lines[0].chr)) {
            if(debug)
                cerr << "heyyyy77"<<lines[0].pos<<" "<<lastchecked<<cytosines[cytosines.size()-1].chr<<endl;
            lastchecked = lines[0].pos;
        }
    }
    if(debug)
        cerr << "lastcheck  :"<<lastchecked<<"  "<<lines[0].pos<<endl;
    if(lastchecked == 341471)
        cerr <<lines[0].seq_string<<endl;
    long refPos = chrom[ChromIndex(lines[0].chr)].chrStart;
    if(debug)
        cerr << "size :"<<lines[0].seq_string<<"   "<< lines[0].chr<<endl;
    int offset = 1;
    if (lastchecked == lines[0].pos) {
        offset = 0;
    }

    for(long i = lastchecked + offset; i< (lines[0].seq_string.size()+lines[0].pos) ;i++){
        if(toupper(reference[refPos + i]) == 'C'){
            //cout<<refPos + i<<"     ref:"<<reference[refPos + i]<<"i:   "<<i<<endl;
            //cout << "rrrr pos:"<<lines[0].pos <<"    at i:"<<lines[0].seq_string[i - lines[0].pos+1] <<"   "<<lines[0].chr<<endl;
            setPointer(i , lines[0].chr);
            cytosines.push_back(Cytosine(i ,lines[0].chr));
            if(debug)
                cerr<<"scn "<<secondPointer<<endl;


            for(int j=0 ; j< secondPointer ; j++){
                //cout<<"seq_string[ i ].pos]  "<<lines[j].seq_string[i - lines[j].pos+1]<<endl;
                if(!(lines[j].pos < (lines[0].pos+lines[0].seq_string.size())))
                   break;
                int relative_pos = lines[j].pos + lines[j].seq_string.size();
                if(i < relative_pos){
                    if (lines[j].seq_string[i - lines[j].pos+1] == 'C') {
                        cytosines[cytosines.size()-1].methylated++;
                    }
                    else if(lines[j].seq_string[1 + i - lines[j].pos] == 'T'){
                        cytosines[cytosines.size()-1].unmethylated++;
                    }
                }
            }
            int s = cytosines.size()-1;
           // fprintf(stdout, "%s\t%ld\t%d\n",cytosines[s].chr, cytosines[s].pos, cytosines[s].methylated);

        }
    }
    //cout<<"sizeeeeee:       "<<lines.size()<<endl;
    lines.erase(lines.begin());
    if(debug)
    cerr<<"sizeeeeee:       "<<lines.size()<<endl;
    computeMethylation();
}

void convertRead(Line &line){
    //cout <<"before :  "<<line.seq_string<<endl;
    int c=0;
    int j = 0;
    for(int i=0;i< line.seq_string.size();i++){
        if (line.cigar2[i] == 'i') {
            line.seq_string.erase(j,1);
            //cout <<"after:  "<<line.seq_string<<endl;
        }
        else if (line.cigar2[i] == 'd') {
            line.seq_string.insert(j,"M");
            j++;
        }
        else
            j++;
        
    }
    //cerr <<"after convert:  "<<line.seq_string<<endl;
}

void reverseRead(Line &line){////
    if(debug)
        cerr <<"before :  "<<line.seq_string<<endl;
    int i,j=0 ;
    char * copy = new char[line.seq_string.size()];
    //    strcpy(copy , line.seq_string);
    line.seq_string.copy(copy, line.seq_string.size(),0);
    for (i = line.seq_string.size()-1; i >= 0 ; i--){
        if(debug)
        cerr << copy[i]<<"c ";
        if(copy[i] == 'A'){
            line.seq_string[j]='T';
        }
        else if(copy[i] == 'T'){
            line.seq_string[j]='A';
        }
        else if(copy[i] == 'C'){
            line.seq_string[j]='G';
        }
        else if(copy[i] == 'G'){
            line.seq_string[j]='C';
        }
        else if(copy[i] == 'M')
            line.seq_string[j] = 'M';
        if(debug)
            cerr << line.seq_string[j]<<"r  ";
        j++;
    }
    if(debug)
        cerr <<"after  reverse:  \n"<<line.seq_string<<endl;
}


char * convertCigar(char * cigar){
    char * cigar2 = new char[500];
    int pos = 0;
    int value = 0, j , index =0;;
    long read_index = 0;
    char alignType;
    //cerr << "cigar in convert"<<cigar<<endl;
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
    //cerr << "cigar2"<<cigar2<<endl;
    return cigar2;
}

int readSamFile(FILE * samFile){
    int index = 0;
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
    
    int count_to=0;
    int stop = 0;
    while (count_to< READ_SAM_BY) {
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
        
        sscanf(line,"%s\t%d\t%s\t%lld\t%u\t%s\t%s\t%s\t%lld\t%s\t%s\n",qname, &flag, temp.chr, &(temp.pos) ,&mapq, temp.cigar,rnext,pnext, &tlen,seq_string,quality_string);
        string str(seq_string);
        temp.seq_string = str;
        
        char * cigar2 = new char[500];
        //cout<<line<<endl;
        
        if(flag != 4){
            //cerr << "hey11111111  "<<endl;
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
       
        int result = checkGAorCT2(temp);
        if(debug)
        cerr << result<<endl;
        if (result == -1) {
            continue;//////////////////////////////////////////
        }
        if(result == 4)
            temp.strand = '+';
        else{
            temp.strand = '-';
            reverseRead(temp);/////////////////////////////////////////////////////////////////
        }
        count_to++;        
        index++;
        if (count_to==READ_SAM_BY) {
            count_to=0;
        }
        lines.push_back(temp);        
    }//while
    if(stop == -1)
        return -1;
    return 1;
    
}

void readSam_Compute(){
    
    int result = readSamFile(samFile);
    if(debug)
        cerr<<result <<"llll        "<<lines.size()<<endl;
    while (result == 1) {
        computeMethylation();
        result = readSamFile(samFile);
        if(result == -1)
            result = 2;
    }
    if(result != 2)
        computeMethylation();
}


void setPointer(int pos, char * chr ){    
    while(secondPointer < lines.size() && lines[secondPointer].pos <= pos && !strcmp(lines[secondPointer].chr,chr)){
        if (lines[secondPointer].seq_string.empty()) {
            break;
        }
        secondPointer++;
    }
    if (secondPointer == lines.size() && lines[secondPointer].pos <= pos && !strcmp(lines[secondPointer].chr,chr)) {
        //int temp = i;
        int res = readSamFile(samFile);
        if(res != -1)
        setPointer(pos, chr);
        
    }
    else if(secondPointer == lines.size())
        secondPointer--;
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
        if(debug)
        for(int k = 0 ; k<200;k++)
            cerr<<reference[k]<<" ";
        samFile = fopen( argv[2], "r" );
        if ( samFile == 0 )
            printf( "Could not open file\n" );
        
        debug = atoi(argv[3]);
        readSam_Compute();
//        if(debug)
//        for(int k=0;k<cytosines.size();k++)
//            cerr <<"cytosin  "<< cytosines[k].pos<<"   "<<cytosines[k].methylated<<"   "<<cytosines[k].chr<<" "<<cytosines[k].unmethylated<<endl;
//        ofstream myfile;
//        myfile.open ("example.txt");
//        myfile << "Writing this to a file.\n";

        
        for (int i=0; i < cytosines.size(); i++) {
            float methylation_ratio = ((float)cytosines[i].methylated/(cytosines[i].methylated+cytosines[i].unmethylated))*100.0 ;
            if((cytosines[i].methylated+cytosines[i].unmethylated)==0)
                methylation_ratio = 0;
            fprintf(stdout, "%s\t%ld\t%d\t%3f\n",cytosines[i].chr, cytosines[i].pos, cytosines[i].methylated,methylation_ratio);
        }
               // myfile.close();
    }
    
}//main


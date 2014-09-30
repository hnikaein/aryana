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
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include "calculate_methylation_acuracy.h"

using namespace std;

long corrects = 0, incorrects = 0;
float max_diffrence = 0.1;

int main(int argc, char *argv[]) {
    char *referenceName, *annotationFile;
    
    if ( argc < 3 ){
        /* We print argv[0] assuming it is the program name */
        printf( "usage: %s original_methylation_filename computed_methylation_filename difference[optional]\n", argv[0] );
    }
    else
    {
		if(argc >= 4)
			max_diffrence = atof(argv[3]);
        calculate(argv[1], argv[2]);
		cout << "Number of correct ratios: " << corrects << endl;
		cout << "Number of incorrect ratios: " << incorrects << endl;
		cout << "Total: " << corrects+incorrects << endl;
    }
    
}//main

void calculate(char* original, char* computed){
    FILE *orig_fp, *comp_fp;
    orig_fp = fopen(original, "r");
	comp_fp = fopen(computed, "r");
	char chrom1[20], chrom2[20];
	long position1, methylated, position2, count = 0;
	float ratio1, ratio2;
	while (! feof(comp_fp)) {
		// comp_fp >> chrom1 >> position1 >> methylated >> ratio1;
		int n = fscanf(comp_fp, "%s\t%ld\t%ld\t%f", chrom1, &position1, &methylated, &ratio1);
		if (n == EOF) break;
		int no_match = 0;
		while(! feof(orig_fp)){
			// orig_fp >> chrom2 >> position2 >> ratio2;
			long file_pos = ftell(orig_fp);
			n = fscanf(orig_fp, "%s\t%ld\t%f", chrom2, &position2, &ratio2);
			if (n == EOF) break;
			if(strcmp(chrom1, chrom2) == 0 && position1 == position2)
				break;
			if(strcmp(chrom1, chrom2) == 0 && position2 > position1){
				no_match = 1;
				fseek(orig_fp, file_pos, SEEK_SET);
				break;
			}
				
		}
		if(no_match){
			if(ratio1 <= max_diffrence)
				corrects++;
			else
				incorrects++;
			continue;
		}
		if(ratio1 >= ratio2 - max_diffrence && ratio1 <= ratio2 + max_diffrence)
			corrects++;
		else
			incorrects++;
		count++;
		if(count % 10 == 0)
			cerr << ".";
	}
}


//
// Created by miladink on 5/11/19.
//
#include <map>
#include <inttypes.h>
#include<iostream>
#include "probnuc.h"
using namespace std;


/*it is a global variable which is read from vcf file,
  it is a map from a position in reference genome to a struct,
  which hold probabilities for the nucleotides at that position,
  the basic schema is like map(1299 -> [0.1, 0.3, 0.4, 0.2], etc).
  In other words, pos_prob_nuc[1299] = {[0.1, 0.3, 0.4, 0.2]} in which
  the second is of probnuc structure
 */
std::map<uint64_t, probnuc> pos_prob_nuc;
std::map<char, int> nuc_2_int = {{'A', 0},{'C', 1},{'G', 2},{'T', 3}};
std::map<int, char> int_2_nuc = {{0, 'A'},{1, 'C'},{2, 'G'},{3, 'T'}};

char sample_from_probnuc(probnuc pn){
    float cum_prob[4];
    cum_prob[0] = pn.prob[0];
    for(int i = 1; i<4; i++){
        cum_prob[i] = pn.prob[i] + cum_prob[i-1];
    }
    double rand_number = (double) rand() / RAND_MAX;
    for(int i = 0; i<4; i++){
        if (rand_number < cum_prob[i])
            return int_2_nuc[i];
    }
    cerr << "Should not reach here" << endl;
}


char give_probnuc_least_chance(probnuc pn){
    int min_chance_nucleotid = 0;
    for(int i = 1; i< 4; i++){
        if(pn.prob[i] < pn.prob[min_chance_nucleotid])
            min_chance_nucleotid = i;
    }
    return int_2_nuc[min_chance_nucleotid];
}


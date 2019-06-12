//
// Created by miladink on 5/11/19.
#include<map>
#ifndef ARYANA_PROBNUC_H
#define ARYANA_PROBNUC_H
/*this structure is for holding probabilities for the existence of one type of nucleotides
 * in the given position. I mean for example if in position 19 of the reference genome,
 * we have probability 0.2 for A, 0.4 for C, and 0.4 for T, and 0 for G, then in this structure
 * we will store [0.2, 0.4, 0.0, 0.4]. The mentioned probabilities are stored in ACGT order.
 * */
struct probnuc{
    float prob[4];
};

extern std::map<uint64_t, probnuc> pos_prob_nuc;
extern std::map<char, int> nuc_2_int;
extern std::map<int, char> int_2_nuc;

char sample_from_probnuc(probnuc pn);

#endif //ARYANA_PROBNUC_H


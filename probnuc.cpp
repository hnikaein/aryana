//
// Created by miladink on 5/11/19.
//
#include <map>
#include <inttypes.h>
#include "probnuc.h"


/*it is a global variable which is read from vcf file,
  it is a map from a position in reference genome to a struct,
  which hold probabilities for the nucleotides at that position,
  the basic schema is like map(1299 -> [0.1, 0.3, 0.4, 0.2], etc).
  In other words, pos_prob_nuc[1299] = {[0.1, 0.3, 0.4, 0.2]} in which
  the second is of probnuc structure
 */
std::map<uint64_t, probnuc> pos_prob_nuc;
std::map<char, int> nuc_2_int = {{'A', 0},{'C', 1},{'G', 2},{'T', 3}};

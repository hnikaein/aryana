//
// Created by miladink on 6/22/19.
//

#ifndef ARYANA_STR_LIB_H
#define ARYANA_STR_LIB_H

#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

vector <string> vector_tokenize_string(string s, string delimiter);
map <string, string> digest_info(string info);
void remove_parenthesis(string& s);

#endif //ARYANA_STR_LIB_H

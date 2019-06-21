//
// Created by miladink on 6/22/19.
//

#include "str_lib.h"

vector <string> vector_tokenize_string(string s, string delimiter) {
    vector <string> ans;
    size_t pos = 0;
    string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        ans.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    ans.push_back(s); //push the last remained part after delimiter
    return ans;
}

map <string, string> digest_info(string info) {
    map <string, string> ans;

    vector<string> pairs = vector_tokenize_string(info,";");
    for (auto pair : pairs) {
        vector<string> keyval = vector_tokenize_string(pair,"=");
        string key = keyval[0];
        string val="";
        if (keyval.size()>1)
            val = keyval[1];
        ans[key] = val;
    }
    return ans;
}

void remove_parenthesis(string& s){
    s.erase(std::remove(s.begin(), s.end(), '['), s.end());
    s.erase(std::remove(s.begin(), s.end(), ']'), s.end());
}

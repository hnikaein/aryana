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
#include <sstream>
#include "const.h"
#include "probnuc.h"

extern std::map<char, int> nuc_2_int;

using namespace std;

vector<long long> chromPos, chromLen;
map<string, int> chromIndex;
map<int, string> chromName;
map<long long, vector<int>> pos_nuc_freq;
int chromNum = 0;
long long genome_size;

//TODO: write a function to convert string to long long it from its string format
//TODO: change string to char
probnuc calculate_probnuc_from_freqs(long long pos) {
    long long total_freq = 0;
    for (int i = 0; i < 4; i++)
        total_freq += pos_nuc_freq[pos][i];

    probnuc pn = {{0.0, 0.0, 0.0, 0.0}};

    for (int i = 0; i < 4; i++)
        pn.prob[i] = (float) pos_nuc_freq[pos][i] / total_freq;
    return pn;
}

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


bool is_file_readable(const string &toTestFile) {
    ifstream ifile(toTestFile.c_str());
    ifile.seekg(0, std::ios_base::end);//seek to end
    //now get current position as length of file
    long long size = ifile.tellg();
    ifile.close();
    if (size <= 0) {
        return false;
    }
    return true;
}

void ReadGenome(const string &genomeFile) {
    genome_size = 0;
    chromNum = 0;

    if (is_file_readable(genomeFile) == false) {
        cout << "Error reading the genome file: " << genomeFile << endl;
        exit(1);
    }

    cerr << "Reading genome..." << endl;
    char fLineMain[10000];
    FILE *f = fopen(genomeFile.c_str(), "r");
    if (!f) {
        cerr << "Error: Genome file not found or could not be opened" << endl;
        exit(1);
    }
    // read genome file piece by piece, line by line
    while (!feof(f)) {
        if (!fgets(fLineMain, sizeof(fLineMain), f)) break;
        int n = static_cast<int>(strlen(fLineMain)), start = 0;
        while (n > 0 && fLineMain[n - 1] <= ' ') n--;
        fLineMain[n] = 0;
        while (start < n && fLineMain[start] <= ' ') start++;
        if (start >= n) continue;
        char *fLine = fLineMain + start;
        n -= start;

        if (fLine[0] == '>') {
            chromPos.push_back(genome_size);
            if (chromNum > 0) {
                chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
            }

            string name = fLine;
            if (name.find(' ') != string::npos) name = name.substr(1, name.find(' ') - 1);
            else name = name.substr(1, name.size() - 1);
            chromIndex[name] = chromNum;
            chromName[chromNum] = name;
            chromNum++;
        } else {
            genome_size += n;
        }
    }
    chromPos.push_back(genome_size);
    chromLen.push_back(chromPos[chromNum] - chromPos[chromNum - 1]);
    fclose(f);
}


long long get_global_pose_from_chromNum_and_index(string chromNumStr, long long local_position) {
    int chromIdx = chromIndex.at(chromNumStr);
    long long global_position = 0;
    for (int i = 0; i < chromIdx; i++) {
        global_position += chromLen[i];
    }
    global_position += local_position;

    return global_position;
}

void print_nuc_probs(ofstream &fout, long long &index_in_aryana, probnuc &pn) {
    fout << index_in_aryana << " ";
    for (int i = 0; i < 3; i++) {
        fout << pn.prob[i] << " ";
    }
    fout << pn.prob[3] << endl;//fout last line separately, so we do not leave extra space at the end of sentence
}

void update_pos_nuc_freq(long long pos, string ref, string alt, long long ref_freq, long long alt_freq) {
    int ref_int = nuc_2_int[ref[0]];
    int alt_int = nuc_2_int[alt[0]];
    pos_nuc_freq[pos].resize(4);
    pos_nuc_freq[pos][ref_int] += ref_freq;
    pos_nuc_freq[pos][alt_int] += alt_freq;
}

int vcf_simplifier(int argc, char *argv[]) {
    string genome_file_path = string(argv[1]);
    string vcf_file_path = string(argv[2]);
    string output_file_path = vcf_file_path + ".simple";

    cout << "vcf file path is " << vcf_file_path << endl;
    cout << "genome file path is " << genome_file_path << endl;

    // process genomeFile to compute chromPos, chromLen, chromIndex, etc.
    ReadGenome(genome_file_path);

    ifstream fin;
    fin.open(vcf_file_path.c_str());
    ofstream fout;
    fout.open(output_file_path.c_str());

    if (!fin) {
        cout << "Cannot open provided vcf file :" << vcf_file_path << endl;
        return 0;
    }

    /*read vcf file line by line, simplifying it and writing the result to a file*/

    string line;
    int line_counter = 0;
    while (getline(fin, line)) {
        if (line == "" || line[0] == '#')
            continue;

        line_counter++;
        istringstream iss(line);
        string chromNum, id, ref, alt, quality, filter, info;
        long long pos;

        iss >> chromNum >> pos >> id >> ref >> alt >> quality >> filter >> info;


        long long global_index_in_genome = get_global_pose_from_chromNum_and_index(chromNum, pos);

        map<string, string> info_map = digest_info(info);
        //TODO just use if AF is not provided?
        long long ref_freq = stoi(info_map["AC"]);
        long long alt_freq = stoi(info_map["AN"]);
        update_pos_nuc_freq(global_index_in_genome, ref, alt, ref_freq, alt_freq);

    }

    //iterating over accumulated frequencies, calculating probabilities
    for (auto const& global_pos_freqs : pos_nuc_freq) {
        long long global_pos = global_pos_freqs.first;
        vector<int> freqs = global_pos_freqs.second;


        probnuc pn = calculate_probnuc_from_freqs(global_pos);

        //calculating where in the aryana index the vcf file is pointing at, genome position and reverse(2*genome_length-genome_position)
        //remember that aryana puts chromosome positions one after another
        long long global_revert_pos = 2 * genome_size - global_pos;
        print_nuc_probs(fout, global_pos, pn);
        print_nuc_probs(fout, global_revert_pos, pn);
    }

    cout << "total vcf lines proceeded: " << line_counter << endl;

    fin.close();
    return 0;
}
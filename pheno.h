#ifndef PHENO_H
#define PHENO_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include "bedFile.h"
#include "htslib/sam.h"
#include <unordered_map>
#include <unordered_set>

using namespace std;

double computeMedian(vector<int>& vec);
double calculateMean(const vector<double>& vec);
double calculateStdDev(const vector<double>& vec);
int pheno_output(int argc,char* argv[]);

bool haveCommonElements(const std::set<string>& set1, const std::set<string>& set2);
map<string,vector<string>> cluster(set<string> intron_set);
void print(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map,string out_file);
void print_total(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map, string refined_output);
pair<map<int,set<string>>,map<int,set<string>>> site_anno(set<int> site, map<int,set<int>> site_corres);
void site_detection(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map, vector<string> sample_name, string inclu_exclu_out, string site_readcount);
int pheno_group(int argc,char* argv[]);



#endif
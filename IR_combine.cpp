#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<vector>
#include<map>
#include<numeric>
#include<utility>
#include "pheno.h"
#include <iomanip>

using namespace std;

map<string,vector<int>> obtain_splice_intron(vector<string> sample,vector<string> site,string intron_file){
    ifstream fin;
    fin.open(intron_file);
    string line,site_name,site_other,clu_other,clu,chr,chr_other,strand,strand_other,start,end,donor_site,acceptor_site;
    map<string,vector<int>> site_used; 
    bool donor_contain,acceptor_contain;
    while(getline(fin,line)){
        clu=line.substr(0,line.find(' '));
        clu_other=line.substr(line.find(' ')+1);
        site_other=clu_other.substr(clu_other.find(' ')+1);
        site_name=clu_other.substr(0,clu_other.find(' '));
        string num_site,tmp;
        chr=site_name.substr(0,site_name.find(':'));
        chr_other=site_name.substr(site_name.find(':')+1);
        strand=chr_other.substr(0,chr_other.find(':'));
        strand_other=chr_other.substr(chr_other.find(":")+1);
        start=strand_other.substr(0,strand_other.find(':'));
        end=strand_other.substr(strand_other.find(':')+1);
        donor_site=chr+":"+strand+":"+start;
        acceptor_site=chr+":"+strand+":"+end;
        if(site_used.find(donor_site)==site_used.end()){
            site_used[donor_site]=vector<int>(sample.size(),0);
        }
        if(site_used.find(acceptor_site)==site_used.end()){
            site_used[acceptor_site]=vector<int>(sample.size(),0);
        }
        for(int i =0;i<sample.size()-1;i++){
            num_site=site_other.substr(0,site_other.find(' '));
            tmp=site_other.substr(site_other.find(' ')+1);
            site_other=tmp;
            site_used[donor_site][i]=site_used[donor_site][i]+stoi(num_site);
            site_used[acceptor_site][i]=site_used[acceptor_site][i]+stoi(num_site);
        }
        site_used[donor_site][sample.size()-1]=site_used[donor_site][sample.size()-1]+stoi(site_other);
        site_used[acceptor_site][sample.size()-1]=site_used[acceptor_site][sample.size()-1]+stoi(site_other);
    }
    fin.close();
    return site_used;
}

map<string,vector<int>> read_in_file(vector<string> sample,vector<string> site,string file_pos,map<string,vector<int>> site_used){
    string line;
    map<string,vector<int>> site_unused=site_used;
    
    for (auto& [key, vec] : site_unused) {
    for (auto& elem : vec) {
        elem = 0;}
    }

    for(int i=0;i<sample.size();i++){
        ifstream fin;
        fin.open(file_pos+"/"+sample[i]+".nonsplit");
        while(getline(fin,line)){
            string chr=line.substr(0,line.find(":"));
            string tmp1=line.substr(line.find(":")+1);
            string pos=tmp1.substr(0,tmp1.find(":"));
            string tmp2=tmp1.substr(tmp1.find(":")+1);
            string strand=tmp2.substr(0,tmp2.find(":"));
            string site_name=chr+":"+strand+":"+pos;
            site_unused[site_name][i]+=1;
        }
        fin.close();
    }
    return site_unused;
}


map<string,vector<string>> filter_output(map<string,vector<int>> site_used,map<string,vector<int>> site_unused){
    for (auto it = site_unused.begin(); it != site_unused.end(); ) {
    int sum = 0;
    for (const auto& s : it->second) {
        sum += s;   // convert string to int
    }

    if (sum == 0) {
        it = site_unused.erase(it); // erase returns next iterator
    } else {
        ++it;
    }

}
    map<string,vector<string>> site_combined;
    for (auto& [key, vec] : site_unused) {
    for(int i=0;i<vec.size();i++){
        string site_used_num=to_string(site_used[key][i]);
        string site_total_num=to_string(site_used[key][i]+site_unused[key][i]);
        string combined_num=site_used_num+":"+site_total_num;
        site_combined[key].push_back(combined_num);
    }
    }
    return site_combined;
}

void print_IR(map<string,vector<string>> output_site,string out_file_prefix,vector<string> sample){
    ofstream fout(out_file_prefix+"_IR.site");
    for(int i =0;i<sample.size()-1;i++){
        fout<<sample[i]<<" ";
    }
    fout<<sample[sample.size()-1]<<endl;
    for(const auto & [value,valuevec]: output_site){
        fout<<value<<" ";
        for(int i=0;i<output_site[value].size()-1;i++){
            fout<<output_site[value][i]<<" ";
        }
        fout<<output_site[value][output_site[value].size()-1]<<endl;
    }
    fout.close();
}

int IR_combine(int argc,char* argv[]){
    string sample_file, file_pos, site_list, intron_file, out_file_prefix;

// Parse options
    int c;
    while ((c = getopt(argc, argv, "hs:f:l:i:o:")) != -1) {
    switch (c) {
        case 'h':
            cout << "Usage: " << argv[0] << " [options]\n"
                 << "Description: Combine single intron cluster site (non-split read combine)\n\n"
                 << "Options:\n"
                 << "  -h          Show this help message and exit\n"
                 << "  -s <file>   Sample file (list of sample names)\n"
                 << "  -f <file>   Non-split read position deposit file\n"
                 << "  -l <file>   Site list file\n"
                 << "  -i <file>   Single intron cluster file\n"
                 << "  -o <prefix> Output file prefix\n";
            exit(0);
        case 's':
            sample_file = string(optarg);
            break;
        case 'f':
            file_pos = string(optarg);
            break;
        case 'l':
            site_list = string(optarg);
            break;
        case 'i':
            intron_file = string(optarg);
            break;
        case 'o':
            out_file_prefix = string(optarg);
            break;
        case '?':
        default:
            throw runtime_error("Error parsing inputs!\n\n");
    }
}

// Validate required arguments
if (sample_file.empty() || file_pos.empty() || site_list.empty() ||
    intron_file.empty() || out_file_prefix.empty()) {
    cerr << "Error: all options -s, -f, -l, -i, -o are required.\n"
         << "Run with -h for usage.\n";
    return 1;
}
    vector<string> sample;
    
    ifstream fin;
    fin.open(sample_file);
    string line;
    while(getline(fin,line)){
        sample.push_back(line);
    }
    fin.close();
    vector<string> site;
    fin.open(site_list);
    while(getline(fin,line)){
        site.push_back(line);
    }
    fin.close();
    map<string,vector<int>> site_used=obtain_splice_intron(sample,site,intron_file);
    map<string,vector<int>> site_unused=read_in_file(sample,site,file_pos,site_used);
    map<string,vector<string>> output_site=filter_output(site_used,site_unused);
    print_IR(output_site,out_file_prefix,sample);
    return 0;
}

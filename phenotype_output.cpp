#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>
#include<cmath>
#include "pheno.h"
using namespace std;
using Mylist = vector<int>;

double computeMedian(vector<int>& vec) {
    // Step 1: Sort the vector
    sort(vec.begin(), vec.end());

    // Step 2: Find the median
    int size = vec.size();
    if (size % 2 == 0) {
        double num = (vec[size / 2 - 1] + vec[size / 2]) / 2.0;
        int intNum = static_cast<int>(round(num));
        // Even number of elements: take the average of the two middle elements
        return intNum;
    } else {
        // Odd number of elements: take the middle element
        return vec[size / 2];
    }
}

double calculateMean(const vector<double>& vec) {
    double sum = 0.0;
    for (double num : vec) {
        sum += num;
    }
    return sum / vec.size();
}

// Function to calculate the standard deviation
double calculateStdDev(const vector<double>& vec) {
    double mean = calculateMean(vec);
    double variance = 0.0;

    for (double num : vec) {
        variance += (num - mean) * (num - mean);
    }
    variance /= vec.size();

    return sqrt(variance);
}

int pheno_output(int argc,char* argv[]){
    string read_count_file, out_file, Prop_file;
    double sd = 0.0, na_prop = 0.0;
    int c;
    while ((c = getopt(argc, argv, "hr:o:p:s:n:")) != -1) {
    switch (c) {
        case 'h':
            cout << "Usage: " << argv[0] << " [options]\n"
                 << "Description: Filter sites with high sparsity and low variance\n\n"
                 << "Options:\n"
                 << "  -h          Show this help message and exit\n"
                 << "  -r <file>   Read count file\n"
                 << "  -o <file>   Output file\n"
                 << "  -p <file>   Proportion file\n"
                 << "  -s <float>  Standard deviation threshold\n"
                 << "  -n <float>  NA proportion threshold\n";
            exit(0);
        case 'r':
            read_count_file = string(optarg);
            break;
        case 'o':
            out_file = string(optarg);
            break;
        case 'p':
            Prop_file = string(optarg);
            break;
        case 's':
            sd = stod(string(optarg));
            break;
        case 'n':
            na_prop = stod(string(optarg));
            break;
        case '?':
        default:
            throw runtime_error("Error parsing inputs!\n\n");
    }}
    if (read_count_file.empty() || out_file.empty() || Prop_file.empty()) {
    cerr << "Error: all options -r, -o, -p, -s, -n are required.\n"
         << "Run with -h for usage.\n";
    return 1;}
   
    ifstream fin;
    fin.open(read_count_file);
    string line,sample,tmp_line;
    string splice_site,readcount;
    string chr_tmp;
    int sp_count,to_count;
    int i=-1;
    vector<string> splice_name, splice_site_;
    map<string,vector<int>> splice_;
    map<string,vector<int>> total_;
    while(getline(fin,line)){
        if(i==-1){
             while(line.find(' ')<100000000){
                sample=line.substr(0,line.find(' '));
                splice_name.push_back(sample);
                tmp_line=line.substr(line.find(' ')+1);
                line=tmp_line;
            }
            sample=line.substr(0,line.find('\n'));
            splice_name.push_back(sample);
        }
        else{
            splice_site=line.substr(0,line.find(' '));
            tmp_line=line.substr(line.find(' ')+1);
            line=tmp_line;
            chr_tmp=splice_site.substr(0,splice_site.find(':'));
            splice_[splice_site]=Mylist();
            total_[splice_site]=Mylist();
            while(line.find(' ')<100000000){
                readcount=line.substr(0,line.find(' '));
                sp_count=stoi(readcount.substr(0,readcount.find(':')));
                to_count=stoi(readcount.substr(readcount.find(':')+1));
                splice_[splice_site].push_back(sp_count);
                total_[splice_site].push_back(to_count);
                tmp_line=line.substr(line.find(' ')+1);
                line=tmp_line;
            }
                readcount=line.substr(0,line.find('\n'));
                sp_count=stoi(readcount.substr(0,readcount.find(':')));
                to_count=stoi(readcount.substr(readcount.find(':')+1));
                splice_[splice_site].push_back(sp_count);
                total_[splice_site].push_back(to_count);
                splice_site_.push_back(splice_site);
            }
            i++;
        }
    cout<<i<<" splice sites read in"<<endl;
   
    ofstream fout(out_file);
    ofstream fout1(Prop_file);
    string key;
    int total;
    for(int i=0;i<splice_name.size();i++){
        fout<<splice_name[i]<<" ";
        fout1<<splice_name[i]<<" ";
    }
    fout<<endl;
    fout1<<endl;
    for(map<string,vector<int>>::iterator it=splice_.begin();it!=splice_.end();++it){
        key=it->first;
        cout<<key<<endl;
        vector<int> na_pos;
        vector<int> non_na;
        vector<int> pre_value;
        vector<int> post_value;
        vector<double> pre_post;
        for(int s=0;s<splice_[key].size();s++){
            total=total_[key][s];
            if(total==0){na_pos.push_back(s);}
            else{
                 non_na.push_back(s);
                 pre_value.push_back(splice_[key][s]);
                 post_value.push_back(total_[key][s]);
                 double pre=splice_[key][s];
                 double post=total;
                 pre_post.push_back(pre/post);
            }
        }
        int median_pre,median_post;
        double size1=na_pos.size();
        double size2=splice_[key].size();
        if(size1/size2<na_prop){
            cout<<na_pos.size()<<"\t"<<splice_[key].size()<<"\t"<<size1/size2<<endl;
            median_pre=computeMedian(pre_value);
            median_post=computeMedian(post_value);
            cout<<median_pre<<"\t"<<median_post<<endl;
            double standard_deviation=calculateStdDev(pre_post);
            cout<<standard_deviation<<endl;
            if(standard_deviation>sd){
                fout<<key<<" ";
                fout1<<key<<" ";
                for(int s=0;s<splice_[key].size();s++){
                    total=total_[key][s];
                    if(total>0){
                        double pre=splice_[key][s];
                        double post=total;
                        if(s<splice_[key].size()-1){
                            fout<<splice_[key][s]<<":"<<to_string(total)<<" ";
                            fout1<<to_string(pre/post)<<" ";}
                        else{
                            fout<<splice_[key][s]<<":"<<to_string(total);
                            fout1<<to_string(pre/post);
                        }
                    }
                    else{
                        double pre=median_pre;
                        double post=median_post;
                        if(s<splice_[key].size()-1){
                            fout<<to_string(median_pre)<<":"<<to_string(median_post)<<" ";
                            fout1<<to_string(pre/post)<<" ";}
                        else{
                            fout<<to_string(median_pre)<<":"<<to_string(median_post);
                            fout1<<to_string(pre/post);
                        }
                    }
                }
                fout<<endl;
                fout1<<endl;
            }
        }
    }
    return 0;
}

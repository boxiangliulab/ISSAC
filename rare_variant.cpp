#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <random>
#include <time.h>
#include <ctime>
#include <string>
#include <cstdlib>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include <unordered_map>
#include <unordered_set>
#include <getopt.h>
#include "eigen_decom.h"
#include "nlopt.hpp"
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <cmath>
#include <gsl/gsl_cdf.h>
#include <omp.h>

using namespace std;

vector<string> split(const string& s, char delimiter) {
    vector<string> tokens;
    string token;
    stringstream ss(s);
    while (getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

int rare_mapping_parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    string vcf_, chr_, output_file_,site_list_,common_sample_file,path_,variant_id,window_,output_path_, threshold_, X_, annot_;
    while((c = getopt(argc, argv, "hs:o:c:v:x:w:p:m:t:a:")) != -1) {
        switch(c) {
            case 'h':
                cout<<"help"<<endl;
            case 's':
                site_list_ = string(optarg);
                break;
            case 'm':
                common_sample_file = string(optarg);
                break;
            case 'c':
                chr_ = string(optarg);
                break;
            case 'v':
                vcf_ = string(optarg);
                break;
            case 'x':
                X_ = string(optarg);
                break;
            case 'w':
                window_ = string(optarg);
                break;
            case 'p':
                path_ = string(optarg);
                break;
            case 'o':
                output_path_ = string(optarg);
                break;
            case 't':
                threshold_ = string(optarg);
                break;
            case 'a':
                annot_ = string(optarg);
                break;
            case '?':
            default:
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    cerr << "Site list: " << site_list_ << endl;
    cerr << "Path: " << path_ << endl;
    cerr << "Output path: " << output_path_ << endl;
    cerr << "chromosome: " << chr_ << endl;
    cerr << "vcf: " << vcf_ << endl;
    cerr << "covariate: " << X_ << endl;
    cerr << "annotation: " << annot_ << endl;
    cerr << endl;

    //read in splice site
    ifstream fin;
    fin.open(site_list_);
    string line,chrom,tmp_line;
    vector<string> sitelist;
    while(getline(fin,line)){
        chrom=line.substr(0,line.find(":"));
        sitelist.push_back(line);
    }
    fin.close();
    cout<<"site list read in"<<endl;
    //read in common_sample
    fin.open(common_sample_file);
    vector<string> common_sample;
    while(getline(fin,line)){
        common_sample.push_back(line);
    }
    fin.close();
    cout<<"sample read in"<<endl;
    //
    fin.open(variant_id);
    vector<string> variant_id_;
    while(getline(fin,line)){
        variant_id_.push_back(line);
    }
    fin.close();
    cout<<"variant read in"<<endl;

    string tmp_1,tmp_2;
    QTL_mapping it;
    it.set_PC(X_);
    //set PC name
    vector<Eigen::VectorXd> test_PC=it.read_in_PC();
    int row = test_PC[0].size();
    int col = test_PC.size();
    cout<<row<<"\t"<<col<<endl;
    Eigen::MatrixXd PCXd(row,col);

    for(int i = 0;i<col;i++){
        PCXd.col(i) = test_PC[i];
    }
    Eigen::MatrixXd X_cons = it.addConstant(PCXd); 
    string tmp_val;
    double val;
    //extract variant annotation
    vector<string> annot_split = split(annot_, ',');

    string chr_site;
    int genotype_num;
    for(int s=0;s<sitelist.size();s++){
        chr_site=sitelist[s].substr(0,sitelist[s].find(":"));
        if(chr_site==chr_){
        ifstream fin1;
        fin1.open(path_+"/"+sitelist[s]+".middle");
        cout<<sitelist[s]<<endl;
        int line_num=0;
        string site_name, dispersion_string;
        double dispersion;
        vector<Eigen::VectorXd> vectors_val_PC;
        bool if_nan=false;
        string dispersion_post;
        while(getline(fin1,line)){
            if(line_num==0){
                cout<<line<<endl;
                tmp_1 = line.substr(line.find("\t")+1);
                site_name=tmp_1.substr(0,tmp_1.find("\t"));
                cout<<site_name<<endl;
                cout<<tmp_1<<endl;
                dispersion_post=tmp_1.substr(tmp_1.find("\t")+1);
                dispersion_string=dispersion_post.substr(0,dispersion_post.find("\t"));
                if(dispersion_string=="nan"){
                    if_nan=true;
                    break;
                }
                dispersion=stod(dispersion_string);
                cout<<dispersion<<"\t"<<dispersion_string<<endl;
            }
            else{
                Eigen::VectorXd val_PC(0);
                tmp_line=line.substr(line.find('\t')+1);
                line=tmp_line;
                while(line.find('\t')<100000000){
                tmp_val=line.substr(0,line.find('\t'));
                val = stod(tmp_val);
                val_PC.conservativeResize(val_PC.size() + 1);
                val_PC(val_PC.size()-1) = val;
                tmp_line=line.substr(line.find('\t')+1);
                line=tmp_line;
            }
            vectors_val_PC.push_back(val_PC);
            }
            line_num++;
        }
        fin1.close();
        if(if_nan==true)continue;
        cout<<"Finish read in middle file"<<endl;
        //random initialize eta,beta_hat,u_hat,mat;
        Eigen::VectorXd eta, beta_hat, u_hat = Eigen::VectorXd::Zero(2);
        Eigen::SparseMatrix<double> result(2,2);
        Eigen::MatrixXd mat(3,3); //X_cons(3,3)
        REMLOptimizer tmp(eta, X_cons, beta_hat, u_hat, vectors_val_PC[3], vectors_val_PC[2], result, 10, 1e-3,mat,0);
        string chrom=site_name.substr(0,site_name.find(":"));
        tmp_1 = site_name.substr(site_name.find(":")+1);
        string pos=tmp_1.substr(tmp_1.find(":")+1);
        vector<pair<string, int>> positions = {{chrom, stoi(pos)}};
        genotype_num=tmp.read_in_genotype_rare(vcf_, annot_split, positions, common_sample,chr_,stoi(window_));
        if(genotype_num>0){
        string output_file=output_path_+"/"+site_name+".rare_result";
        tmp.collapse_genotype(stoi(window_));
        tmp.pvalue_beta_sd_compute_rare(output_file, site_name, dispersion,stod(threshold_),vectors_val_PC[0],vectors_val_PC[1],vectors_val_PC[2],vectors_val_PC[3]);}  
    }}
    return 0;
}

int rare_mapping(int argc, char *argv[]) {
    int i = rare_mapping_parse_options(argc,argv);
    return 0;
}
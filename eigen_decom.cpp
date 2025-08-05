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
using Mylist = vector<int>;

int QTL_mapping::parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    while((c = getopt(argc, argv, "hs:o:u:p:g:n:")) != -1) {
        switch(c) {
            case 'h':
                cout<<"help"<<endl;
            case 's':
                splice_file_ = string(optarg);
                break;
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'u':
                output_path_ = string(optarg);
                break;
            case 'p':
                PC_ = string(optarg);
                break;
            case 'g':
                GRM_ = string(optarg);
                break;
            /*case 'c':
                chr_ = string(optarg);
                break;
            case 'v':
                vcf_ = string(optarg);
                break;*/
            case 'n':
                GRM_num_ = string(optarg);
                break;
            case '?':
            default:
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    chr_="NA";
    vcf_="NA";
    cerr << "Phenotype: " << splice_file_ << endl;
    cerr << "PC: " << PC_ << endl;
    cerr << "GRM: " << GRM_ << endl;
    cerr << "GRM_num: " << GRM_num_ << endl;
    cerr << "Output file: " << output_file_ << endl;
    //cerr << "chromosome: " << chr_ << endl;
    //cerr << "vcf: " << vcf_ << endl;
    cerr << endl;
    return 0;
}

//read in phenotype
int QTL_mapping::read_in_splice(){
    ifstream fin;
    fin.open(splice_file_);
    if (!fin.is_open()) {
            cerr << "Failed to open the phenotype file." << endl;
            return 1; // Or handle the error as appropriate
            }
    int i=-1;
    string line,sample,tmp_line;
    string splice_site,readcount;
    string chr_tmp;
    int sp_count,to_count;
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
            chr_=chr_tmp; // avoid only one chr could be read in
            if(chr_tmp==chr_){
                splice_[splice_site]=Mylist();
                splice_unsplice_[splice_site]=Mylist();
                while(line.find(' ')<100000000){
                    readcount=line.substr(0,line.find(' '));
                    sp_count=stoi(readcount.substr(0,readcount.find(':')));
                    to_count=stoi(readcount.substr(readcount.find(':')+1));
                    splice_[splice_site].push_back(sp_count);
                    splice_unsplice_[splice_site].push_back(to_count);
                    tmp_line=line.substr(line.find(' ')+1);
                    line=tmp_line;
                }
                splice_site_.push_back(splice_site);
                cout<<"addi"<<line<<endl;
                //readcount=line.substr(0,line.find(' '));
                readcount = line;
                cout<<readcount<<endl;
                sp_count=stoi(readcount.substr(0,readcount.find(':')));
                to_count=stoi(readcount.substr(readcount.find(':')+1));
                //cout<<splice_site<<"\t"<<sp_count<<"\t"<<to_count<<endl;
                splice_[splice_site].push_back(sp_count);
                splice_unsplice_[splice_site].push_back(to_count);
            }
        }
        i++;
    }
    cout<<i<<" splice sites read in"<<endl;
    return 0;
}

//read in PC file
vector<Eigen::VectorXd> QTL_mapping::read_in_PC(){
    ifstream fin;
    fin.open(PC_);
    if (!fin.is_open()) {
            cerr << "Failed to open the PC file." << endl;
             // Or handle the error as appropriate
            }
    string line,tmp_val,sample,tmp_line;
    double val;
    vector<Eigen::VectorXd> vectors_val_PC;
    int i = -1;
    while(getline(fin,line)){
         if(i==-1){
             while(line.find('\t')<100000000){
                sample=line.substr(0,line.find('\t'));
                PC_name.push_back(sample);
                tmp_line=line.substr(line.find('\t')+1);
                line=tmp_line;
            }
            sample=line.substr(0,line.find('\n'));
            PC_name.push_back(sample);
        }
        else{
            Eigen::VectorXd val_PC(0);
            while(line.find('\t')<100000000){
                tmp_val=line.substr(0,line.find('\t'));
                val = stod(tmp_val);
                val_PC.conservativeResize(val_PC.size() + 1);
                val_PC(val_PC.size()-1) = val;
                tmp_line=line.substr(line.find('\t')+1);
                line=tmp_line;
            }
            tmp_val=line.substr(0,line.find('\n'));
            val = stod(tmp_val);
            val_PC.conservativeResize(val_PC.size() + 1);
            val_PC(val_PC.size()-1) = val;
            vectors_val_PC.push_back(val_PC);
        }
        i++;
    }
    PC_mat_ = vectors_val_PC;
    return vectors_val_PC;
}


// Function to generate a random sparse symmetric positive-definite matrix
Eigen::SparseMatrix<double> QTL_mapping::read_in_GRM() {
    int n = stoi(GRM_num_);
    Eigen::SparseMatrix<double> GRM(n,n);
    ifstream fin(GRM_);
    if (!fin.is_open()) {
            cerr << "Failed to open the GRM file." << endl;
             // Or handle the error as appropriate
            }
    string line,sample,tmp_line,tmp_val;
    double val;
    int i=-1;
    //Check the length of GRM;
    while(getline(fin,line)){
        if(i==-1){
            for(int j=0;j<n;j++){
                sample=line.substr(0,line.find(' '));
                GRM_name.push_back(sample);
                if(j<n-1){
                    tmp_line=line.substr(line.find(' ')+1);
                    line=tmp_line;}
            }
        }
        else{
            sample=line.substr(0,line.find(' '));
            tmp_line=line.substr(line.find(' ')+1);
            line=tmp_line;
             for(int j=0;j<n;j++){
                tmp_val = line.substr(0,line.find(' '));
                val = stod(tmp_val);
                if(val>0.6){
                GRM.insert(i,j) = val;}
                if(j<n-1){
                    tmp_line=line.substr(line.find(' ')+1);
                    line=tmp_line;
                }
            }
        }
        i++;
    }
    return GRM;
}

Eigen::SparseMatrix<double> QTL_mapping::identity_matrix(){
    Eigen::SparseMatrix<double> identity(common_name.size(),GRM_name.size());
    for(int i=0;i<common_name.size();i++){
        string tmp_ = common_name[i].substr(0,common_name[i].find(':'));
        for(int j=0;j<GRM_name.size();j++){
            if(tmp_==GRM_name[j]){
                identity.insert(i,j) = 1;
                break;
            }
        }
    }
    return identity;
}

//obtain common name for all input files and use common name for downstream analysis
int QTL_mapping::findElementIndex(vector<string> vec, string value){
    for (int i = 0; i < vec.size(); ++i) {
        if (vec[i] == value) {
            return i; // Return the index if found
        }
    }
    return -1;
}

void QTL_mapping::obtain_common_name(){
    for(int i=0;i<splice_name.size();i++){
        for(int j=0;j<PC_name.size();j++){
            if(splice_name[i]==PC_name[j]){
                common_name.push_back(splice_name[i]);
                break;
            }
        }
    }

    //use common name for downstream analysis
    vector<Eigen::VectorXd> PC_mat_tmp = PC_mat_;
    map<string,vector<int>> splice_tmp = splice_;
    map<string,vector<int>> splice_unsplice_tmp = splice_unsplice_;

    
    for(auto it = splice_tmp.begin();it!=splice_tmp.end();++it){
        string key = it->first;
        splice_tmp[key].clear();
        splice_unsplice_tmp[key].clear();
    }
    cout<<"Finish clear"<<endl;

    for(int i = 0;i<PC_mat_tmp.size();i++){
        PC_mat_tmp[i].resize(0);
    }
    for(int i = 0;i<common_name.size();i++){
        //adjust PC
        int PC_pos = findElementIndex(PC_name,common_name[i]);
        for(int j=0;j<PC_mat_.size();j++){
            PC_mat_tmp[j].conservativeResize(PC_mat_tmp[j].size() + 1);
            PC_mat_tmp[j](PC_mat_tmp[j].size()-1) = PC_mat_[j][PC_pos];
        }
        //adjust splice_name
        int splice_pos = findElementIndex(splice_name,common_name[i]);
        for(const auto& pair: splice_){
        string key = pair.first;
        splice_tmp[key].push_back(splice_[key][splice_pos]);
        splice_unsplice_tmp[key].push_back(splice_unsplice_[key][splice_pos]);}
    }
    PC_mat_ = PC_mat_tmp;
    splice_ = splice_tmp;
    splice_unsplice_ = splice_unsplice_tmp;
    cout<<"Finish change"<<endl;
}

Eigen::VectorXd REMLOptimizer::compute_V_inv_X(Eigen::SparseMatrix<double> result, Eigen::VectorXd b, int n, int max_iteration, double tol) {
    Eigen::VectorXd x(n);
    // Set up the Conjugate Gradient solver with a diagonal preconditioner
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    // Compute the decomposition of A
    cg.compute(result);
    cg.setMaxIterations(max_iteration);
    cg.setTolerance(tol);
    // Solve the system A * x = b using PCG
    
    x = cg.solve(b);

    // Check for convergence
    //if (cg.info() == Eigen::Success) {
    //    std::cout << "Conjugate Gradient converged successfully." << std::endl;
    //} else {
    //    std::cout << "Conjugate Gradient failed to converge." << std::endl;
    //}
    // Output the number of iterations and the solution's first 10 elements
    //std::cout << "Number of iterations: " << cg.iterations() << std::endl;
    //std::cout << "Estimated error: " << cg.error() << std::endl;
    //std::cout << "First 10 elements of the solution vector x: " << x.head(10).transpose() << std::endl;
    return x;
}


double REMLOptimizer::logAbsDeterminant(const Eigen::SparseMatrix<double>& V) {
    // Perform LU decomposition for the sparse matrix
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(V);  // Analyze the sparsity pattern
    solver.factorize(V);       // Factorize the matrix

    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1.0;  // Handle decomposition failure
    }
    double log_det = solver.logAbsDeterminant();
    // Compute log(|det(V)|) as the sum of log of the absolute value of the diagonal elements of U
    return log_det;
}


void QTL_mapping::Bino_GLMM(string site,Eigen::SparseMatrix<double> result){
    // initialize: obtain beta, eta(predicted value)
    vector<int> splice_num = splice_[site];
    vector<int> total_num = splice_unsplice_[site];
    //Transform splice_num, total num to VectorXd
    Eigen::VectorXd splice_numXd(splice_[site].size());
    Eigen::VectorXd total_numXd(splice_unsplice_[site].size());
    // Step 3: Copy the values from std::vector to Eigen::VectorXd
    for (size_t i = 0; i < splice_[site].size(); ++i) {
        splice_numXd(i) = splice_[site][i];
        total_numXd(i) = splice_unsplice_[site][i];
    }

    //Transform PC_mat to MatrixXd
    int row = PC_mat_[0].size();
    int col = PC_mat_.size();
    Eigen::MatrixXd PCXd(row,col);

    for(int i = 0;i<col;i++){
        PCXd.col(i) = PC_mat_[i];
    }

    Eigen::MatrixXd X_cons = addConstant(PCXd);  // X_cons will have an intercept column
    cout<<"X_cons"<<X_cons.cols()<<"\t"<<X_cons.rows()<<endl;
    // Perform IRLS to obtain beta estimates for the Binomial GLM
    Eigen::VectorXd beta_hat = IRLS(X_cons, splice_numXd, total_numXd, 10, 1e-3);
    cout<<"splice site:"<<splice_numXd.size()<<"\t"<<total_numXd.size()<<endl;
    //Eigen::VectorXd beta_hat = Eigen::VectorXd::Zero(X_cons.cols());
    cout<<beta_hat<<endl;
    Eigen::VectorXd eta = X_cons * beta_hat;
    //cout<<eta<<endl;
    // optimize to find the best estimates for tau
    int iter=20;
    double tol = 1e-6;
    Eigen::MatrixXd mat(3,3);
    mat.setRandom();
    //initialize u_hat
    Eigen::VectorXd u_hat = Eigen::VectorXd::Zero(splice_numXd.size());
    REMLOptimizer tmp(eta, X_cons, beta_hat, u_hat, splice_numXd, total_numXd, result, iter, tol,mat,0);
    tmp.update(splice_numXd,total_numXd);
    string chrom = site.substr(0,site.find(":"));
    string tmp_pos = site.substr(site.find(":")+1);
    string pos = tmp_pos.substr(tmp_pos.find(":")+1);
    std::vector<std::pair<std::string, int>> positions = {
        {chrom, stoi(pos)}
    };
    //tmp.read_in_genotype(vcf_, positions, windowsize, common_name);
    tmp.compute(output_path_,site);
}


void QTL_mapping::test(){
    // Generate a large sparse symmetric positive-definite matrix A (3000x3000)
    //read in identity matrix
    Eigen::SparseMatrix<double> A = read_in_GRM();
    cout<<A.rows()<<"\t"<<A.cols()<<endl;
    cout<<"read in GRM success"<<endl;
    
    int max_iteration=100;
    double tol=1e-6;
    windowsize = 1000000;
    //Eigen::VectorXd b(n);
    //b.setRandom();
    //Eigen::VectorXd x = compute_V_inv_X(result, b, n, max_iteration, tol);
    //cout<<x.head(10).transpose()<<endl;
    //double log_det = logAbsDeterminant(result);
    //cout<<log_det<<endl;
    cout<<PC_<<endl;
    vector<Eigen::VectorXd> test_PC=read_in_PC();
    int n = test_PC[0].size();
    int i=read_in_splice();

    obtain_common_name();
    cout<<splice_name.size()<<"\t"<<PC_name.size()<<endl;
    Eigen::SparseMatrix<double> Identity = identity_matrix();

    Eigen::SparseMatrix<double> result = (Identity*(A))*(Identity.transpose());
    //for(int i=0;i<n;i++){
    //    result.coeffRef(i,i) = result.coeffRef(i,i) + 0.3;
    //}
    cout<<"GRM extended"<<endl;
    cout<<result.rows()<<"\t"<<result.cols()<<endl;
    for(int site = 0;site<splice_site_.size();site++){
        cout<<splice_site_[site]<<endl;
        Bino_GLMM(splice_site_[site],result);
    }
}



// Function to add intercept (constant) to the design matrix
Eigen::MatrixXd QTL_mapping::addConstant(const Eigen::MatrixXd& X) {
    int n = X.rows();
    Eigen::MatrixXd X_cons(n, X.cols() + 1);
    X_cons.col(0) = Eigen::VectorXd::Ones(n);  // Intercept column
    X_cons.block(0, 1, n, X.cols()) = X;  // Rest of the columns
    return X_cons;
}

// IRLS function to estimate beta
Eigen::VectorXd QTL_mapping::IRLS(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, const Eigen::VectorXd& total, int maxIter, double tol) {
    int n = X.rows();
    int p = X.cols();
    cout<<X.cols()<<"\t"<<X.rows()<<endl;
    cout<<y.size()<<"\t"<<total.size()<<endl;
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);  // Initialize beta
    Eigen::VectorXd eta, pi,mu, z, XtWy,W;
    Eigen::MatrixXd W_diag, XtW, XtWX;
    Eigen::SparseMatrix<double> XtWX_sparse;
    double pre_beta = 0;
    for (int iter = 0; iter < maxIter; ++iter) {
        // Step 1: Compute eta = X * beta
        eta = X * beta;
        cout<<eta.head(10).transpose()<<endl;
        // Step 2: Compute mu = logit(eta)
        pi = logit(eta);
        mu = pi.array()*total.array();
        // Step 3: Compute the working dependent variable z
        z = eta.array() + (y.array() - mu.array()) / (total.array()* pi.array() * (1 - pi.array()));
        cout<<z.size()<<endl;
        // Step 4: Compute weights (W_diag = mu * (1 - mu))
        W = ((pi.array()*(1-pi.array()))).matrix();
        W_diag = W.asDiagonal();
        //Eigen::MatrixXd W_inv = W_diag.inverse();
        // Convert to sparse matrix
        //Eigen::SparseMatrix<double> W_inv_sparse = W_diag.sparseView();

        // Step 5: Update beta using weighted least squares
        cout<<"pass1"<<endl;
	XtW = X.transpose() * W_diag;
	cout<<XtW.cols()<<"\t"<<XtW.rows()<<endl;
        XtWX = XtW * X;
        XtWX_sparse = XtWX.sparseView();
        XtWy = XtW * z;
        pre_beta = beta[0];
        cout<<"pass2"<<endl;
	cout<<XtWy.head(2)<<endl;
        Eigen::VectorXd x(X.rows());
	cout<<"pass3"<<endl;
    // Set up the Conjugate Gradient solver with a diagonal preconditioner
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    // Compute the decomposition of A
        cg.compute(XtWX_sparse);
        cg.setMaxIterations(maxIter);
        cg.setTolerance(tol);
    // Solve the system A * x = b using PCG
    
        beta = cg.solve(XtWy);  // Solve for new beta
        cout<<"beta\t"<<beta.head(2).transpose()<<endl;
        // Step 6: Check for convergence
        if (abs(beta[0]-pre_beta) < tol) {
            break;  // Convergence reached
        }
        beta = Eigen::VectorXd::Zero(p); // if not converged
    }
    if(isnan(beta[0])){beta = Eigen::VectorXd::Zero(p);}
    return beta;
}


int model(int argc, char *argv[]) {
    QTL_mapping it;
    it.parse_options(argc,argv);
    it.test();
    return 0;
}


int QTL_parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    string vcf_, chr_, output_file_,site_list_,common_sample_file,path_,windowsize,output_path_, threshold_, X_;
    while((c = getopt(argc, argv, "hs:o:c:v:x:p:w:m:t:")) != -1) {
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
            case 'p':
                path_ = string(optarg);
                break;
            case 'w':
                windowsize = string(optarg);
                break;
            case 'o':
                output_path_ = string(optarg);
                break;
            case 't':
                threshold_ = string(optarg);
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
    cerr << endl;

    //read in splice site
    ifstream fin;
    fin.open(site_list_);
    string line,chrom,tmp_line;
    vector<string> sitelist;
    while(getline(fin,line)){
        chrom=line.substr(0,line.find(":"));
        if(chrom==chr_){
            sitelist.push_back(line);
        }
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
    for(int s=0;s<sitelist.size();s++){
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
        tmp.read_in_genotype(vcf_, positions, stoi(windowsize), common_sample);
        string output_file=output_path_+"/"+site_name+".result";
        tmp.pvalue_beta_sd_compute(output_file, site_name, dispersion,stod(threshold_),vectors_val_PC[0],vectors_val_PC[1],vectors_val_PC[2],vectors_val_PC[3]);  
    }
    return 0;
}

int DS_parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    string output_file_,site_list_,common_sample_file,path_,output_path_,output_prefix_,group_file_;
    while((c = getopt(argc, argv, "hs:o:p:m:g:")) != -1) {
        switch(c) {
            case 'h':
                cout<<"help"<<endl;
            case 's':
                site_list_ = string(optarg);
                break;
            case 'm':
                common_sample_file = string(optarg);
                break;
            case 'g':
                group_file_ = string(optarg);
                break;
            case 'p':
                path_ = string(optarg);
                break;
            case 'o':
                output_path_ = string(optarg);
                break;
            case '?':
            default:
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    cerr << "Site list: " << site_list_ << endl;
    cerr << "Path: " << path_ << endl;
    cerr << "Output path: " << output_path_ << endl;
    cerr << "Group file: " << group_file_ << endl;
    cerr << endl;

    //read in splice site
    ifstream fin;
    fin.open(site_list_);
    string line,tmp_line;
    vector<string> sitelist;
    while(getline(fin,line)){
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
    //read in group file
    fin.open(group_file_);
    vector<int> group;
    while(getline(fin,line)){
        group.push_back(stoi(line));
    }
    fin.close();
    //
    string tmp_1,tmp_2;
    string tmp_val;
    double val;
    for(int s=0;s<sitelist.size();s++){
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
        Eigen::MatrixXd mat(3,3),X_cons(3,3);
        REMLOptimizer tmp(eta, X_cons, beta_hat, u_hat, vectors_val_PC[3], vectors_val_PC[2], result, 10, 1e-3,mat,0);
        string output_file=output_path_+"/"+sitelist[s]+".DS_result";
        tmp.DS(output_file, site_name, dispersion,group, vectors_val_PC[0],vectors_val_PC[1],vectors_val_PC[2],vectors_val_PC[3]);
    }
    return 0;
}

int trans_QTL_parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    string vcf_, chr_, output_file_,site_list_,common_sample_file,path_,variant_id,output_path_, threshold_, X_;
    while((c = getopt(argc, argv, "hs:o:c:v:x:w:p:m:t:")) != -1) {
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
                variant_id = string(optarg);
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
    for(int s=0;s<sitelist.size();s++){
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
        tmp.read_in_genotype_trans(vcf_, variant_id_, common_sample,chr_);
        string output_file=output_path_+"/"+site_name+".result";
        tmp.pvalue_beta_sd_compute(output_file, site_name, dispersion,stod(threshold_),vectors_val_PC[0],vectors_val_PC[1],vectors_val_PC[2],vectors_val_PC[3]);  
    }
    return 0;
}


int QTL(int argc, char *argv[]) {
    int i = QTL_parse_options(argc,argv);
    return 0;
}

int DS(int argc, char *argv[]) {
    int i = DS_parse_options(argc,argv);
    return 0;
}


int trans_QTL(int argc, char *argv[]) {
    int i = trans_QTL_parse_options(argc,argv);
    return 0;
}





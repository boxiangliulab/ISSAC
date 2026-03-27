#ifndef EIGEN_DECOM
#define EIGEN_DECOM

#include <algorithm>
#include <iomanip>
#include <iostream>
#include "htslib/sam.h"
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <tuple>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <time.h>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include "nlopt.hpp"
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/hts.h>
#include <gsl/gsl_cdf.h>
#include <getopt.h>
#include <omp.h> 

using namespace std;
using Mylist = vector<int>;
//The class that deals with QTL mapping
class QTL_mapping {
    private:
        //(splice read count) : (unsplice read count + splice read count)
        string splice_file_;
        //File to write output to - optional, write to STDOUT by default
        string output_file_;
        string output_path_;
        //File to PC - optional
        string PC_;
        //File to GRM - optional
        string GRM_;
        //File to vcf.gz
        string vcf_;
        // deposite chr
        string chr_;
        string GRM_num_;
        // iterate times for normalization parameter estimation
        int time_;
        int windowsize;
        // deposite splice read count and total read count
        map<string,vector<int>> splice_;
        map<string,vector<int>> splice_unsplice_;
        vector<Eigen::VectorXd> PC_mat_;
        //deposit splice site name
        vector<string> splice_site_;
        // deposite sample name;
        vector<string> splice_name;
        vector<string> PC_name;
        vector<string> GRM_name;
        vector<string> common_name; // samples used for QTL mapping
    public:
        //Default constructor
        QTL_mapping() {
            splice_file_ = "NA";
            output_file_ = "NA";
            output_path_ = "NA";
            PC_ = "NA";
            GRM_ = "NA";
            chr_ = ".";
            time_=0;
            vcf_ = "NA";
            GRM_num_ = "NA";
        }

        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]); 
        
        //read in GRM file
        Eigen::SparseMatrix<double> read_in_GRM();
        
        //read in PC file
        vector<Eigen::VectorXd> read_in_PC();
        void set_PC(string name){
            PC_ = name;
        }
        //read in splice site phenotype
        int read_in_splice();

        void obtain_common_name();
        // facilitate computation
        // Binomial GLMM model training
        void Bino_GLMM(string site,Eigen::SparseMatrix<double> result,Eigen::SparseMatrix<double> Identity);
        
        Eigen::MatrixXd addConstant(const Eigen::MatrixXd& X);

        Eigen::VectorXd IRLS(const Eigen::MatrixXd& X, const Eigen::VectorXd& y,const Eigen::VectorXd& total, int maxIter, double tol);

        int findElementIndex(vector<string> vec, string value);

        Eigen::VectorXd logit(const Eigen::VectorXd& eta) {
            return 1.0 / (1.0 + (-eta).array().exp());}
        
        Eigen::SparseMatrix<double> identity_matrix();
        //For temporary modify
        void test();

};


class REMLOptimizer {
public:
    // This class wraps your optimization problem

    REMLOptimizer(Eigen::VectorXd& eta, const Eigen::MatrixXd& X, 
                  Eigen::VectorXd& beta_hat,Eigen::VectorXd& u_hat_, Eigen::VectorXd& y_, 
                  Eigen::VectorXd& total_, Eigen::SparseMatrix<double> K, const int max_iteration,
                  const double tol, Eigen::MatrixXd& sigma_inv_X, double tau_) 
        : eta(eta), X(X), beta_hat(beta_hat), u_hat_(u_hat_),y_(y_),total_(total_), K(K), max_iteration(max_iteration), tol(tol), sigma_inv_X(sigma_inv_X), tau_(tau_) {}

    static double reml_criterion(const vector<double> &theta, vector<double> &grad, void *data){
        auto *args = static_cast<std::tuple<Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> *>(data);
        Eigen::VectorXd eta = get<0>(*args);
        Eigen::MatrixXd X = get<1>(*args);
        Eigen::VectorXd beta_hat = get<2>(*args);
        Eigen::SparseMatrix<double> W_inv = get<3>(*args);
        Eigen::SparseMatrix<double> K = get<4>(*args);
        
        double tau2 = std::exp(theta[0]);
        int max_iteration=100;
        double tol=1e-3;
        Eigen::SparseMatrix<double> result = W_inv + tau2 * K;
        Eigen::MatrixXd V_inv_X = Eigen::MatrixXd::Zero(X.rows(),X.cols());
        for(int i=0;i<X.cols();i++){
            Eigen::VectorXd b_i = X.col(i);
            Eigen::VectorXd x(X.rows());
    // Set up the Conjugate Gradient solver with a diagonal preconditioner
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    // Compute the decomposition of A
            cg.compute(result);
            cg.setMaxIterations(max_iteration);
            cg.setTolerance(tol);
    // Solve the system A * x = b using PCG
    
            x = cg.solve(b_i);
            //Eigen::VectorXd x = compute_V_inv_X(result, b_i, X.rows(), max_iteration, tol);
            V_inv_X.col(i) = x;
        }
        //sigma_inv_X = V_inv_X;
        Eigen::SparseMatrix<double> X_t_V_inv_X = (X.transpose() * V_inv_X).sparseView();
        Eigen::VectorXd residual = eta - X * beta_hat;
        Eigen::VectorXd x(result.rows());
    // Set up the Conjugate Gradient solver with a diagonal preconditioner
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::DiagonalPreconditioner<double>> cg;
    // Compute the decomposition of A
        cg.compute(result);
        cg.setMaxIterations(max_iteration);
        cg.setTolerance(tol);
    // Solve the system A * x = b using PCG
    
        Eigen::VectorXd V_inv_residual = cg.solve(residual);
        //Eigen::VectorXd V_inv_residual = compute_V_inv_X(result, residual , result.rows(), max_iteration, tol);
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(result);  // Analyze the sparsity pattern
        solver.factorize(result);       // Factorize the matrix

        if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1.0;  // Handle decomposition failure
        }
        double log_det_V = solver.logAbsDeterminant();
        //double log_det_V = logAbsDeterminant(result);
        //double log_det_X_t_V_inv_X = logAbsDeterminant(X_t_V_inv_X);
        //double log_det_X_t_V_inv_X = 0;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_2;
        solver_2.analyzePattern(X_t_V_inv_X);  // Analyze the sparsity pattern
        solver_2.factorize(X_t_V_inv_X);       // Factorize the matrix

        if (solver_2.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1.0;  // Handle decomposition failure
        }
        double log_det_X_t_V_inv_X = solver_2.logAbsDeterminant();

        double reml_value = -0.5*(log_det_V + log_det_X_t_V_inv_X + residual.transpose()*V_inv_residual);
        cout<<-reml_value<<endl;
        return -reml_value;
    }

    Eigen::VectorXd sigmoid(const Eigen::VectorXd & x){
        return 1.0/(1.0+(-x.array()).exp());
    }

    void update(Eigen::VectorXd& y,Eigen::VectorXd& total){
        
        y_ = y;
        total_ = total;
        double para_init = 0;
        double tau = para_init;
    
        //Eigen::VectorXd beta_hat; // Initialize beta_hat with appropriate values
        Eigen::VectorXd pre_beta = beta_hat;
        double pre_tau = 0;
        double pre_phi = 0;

        //Eigen::MatrixXd newX;  // Initialize newX with your data
    // Example: Eigen::MatrixXd newX = Eigen::MatrixXd::Constant(100, 2, 1.0);

    // Initialize other variables: total, u_hat, y, phi, K, etc.
        Eigen::VectorXd u_hat = Eigen::VectorXd::Zero(y.size());
        //Eigen::MatrixXd K;

    // Initialize pi, mu, eta
        Eigen::VectorXd pi = sigmoid(X * beta_hat + u_hat);
        Eigen::VectorXd mu = pi.array() * total.array();
        //Eigen::VectorXd eta = (X * beta_hat).array() + u_hat.array() + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));

        double phi = 1.0;

        //int iter = 100;  // Set the number of iterations (as in the Python code)
        //double tol = 1e-6; // Tolerance for convergence
        cout<<"Start iteration"<<endl;
        for (int i = 0; i < max_iteration; ++i) {
        // Update pi, mu, eta
        pi = sigmoid(X * beta_hat + u_hat);
        mu = pi.array() * total.array();
        eta = (X * beta_hat).array() + u_hat.array() + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));
        //cout<<eta[0]<<endl; //
        // Compute W and W_inv (inverse of W)
        Eigen::VectorXd W_vec = ((pi.array()*(1-pi.array())*phi)).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();

        Eigen::VectorXd W_inv_vec = (1/(pi.array()*(1-pi.array())*phi)).matrix();
        Eigen::MatrixXd W_inv_diag = W_inv_vec.asDiagonal();
        Eigen::SparseMatrix<double> W_inv_sparse = W_inv_diag.sparseView();

        // Initialize theta (equivalent to np.log([np.var(eta)]))
        double mean_eta = eta.mean();
        double variance_eta = (eta.array()-mean_eta).square().mean();
        double log_var_eta = log(variance_eta);
        cout<<"log_var_eta\t"<<log_var_eta<<endl;
        if(std::isnan(log_var_eta)){
            throw std::invalid_argument("Iteration failed");
        }
        vector<double> theta(1,log_var_eta);
        cout<<"eta"<<"\t"<<eta.head(10).transpose()<<endl;
        cout<<"pi"<<"\t"<<pi.head(10).transpose()<<endl;
        cout<<"y"<<"\t"<<y.head(10).transpose()<<endl;
        cout<<"total"<<"\t"<<total.head(10).transpose()<<endl;
        // Define the optimizer
        //LBFGSpp::LBFGSParam<double> param;
        //param.epsilon = tol;
        //LBFGSpp::LBFGSSolver<double> solver(param);
        //double fx;
        cout<<"Defined successful"<<endl;
        double tau_pre = 0;
        double tau_post = log_var_eta;
        //test using nlopt.hpp
        auto args = std::make_tuple(eta,X, beta_hat,W_inv_sparse,K);
        nlopt::opt opt(nlopt::LN_COBYLA,1);
        std::vector<double> lb(1,-HUGE_VAL);  // Lower bounds for x
        std::vector<double> ub(1,HUGE_VAL);
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);
    // Set the objective function
        opt.set_min_objective(reml_criterion, &args);

    // Set tolerance and stopping criteria (optional)
        opt.set_xtol_rel(1e-6);

        double minf;  // Minimum value of the objective function

        try {
        // Run the optimization
           nlopt::result result = opt.optimize(theta, minf);
           //nlopt::result result = opt.optimize(theta);
           std::cout << "Found minimum at f(x) = " << minf << std::endl;
           std::cout << "Optimal parameters: " << theta[0]  << std::endl;
         } catch (std::exception &e) {
           std::cerr << "NLopt failed: " << e.what() << std::endl;
         }

        //auto rem_optimizer = [&](Eigen::VectorXd& theta, Eigen::VectorXd& grad) -> double {
        //    return reml_criterion(theta, y, eta,X,beta_hat,W_inv_sparse,K); // Implement gradient computation if needed
        //};

        //int niter = solver.minimize(rem_optimizer, theta, fx);*/
        
        // Update tau
        double tau_next = std::exp(theta[0]);

        //deal with abnormal condition

        std::cout << "Iteration " << i << " : " << tau_next << std::endl;
        
        // Convergence check
        pre_tau = tau;
        tau = tau_next;

        // Update sigma and sigma_inv
        Eigen::SparseMatrix<double> sigma = W_inv_sparse + tau * K;
        //Eigen::MatrixXd sigma_inv = sigma.inverse();
        // Update beta_hat
        pre_beta = beta_hat;
        //Eigen::VectorXd sigma_inv_X = compute_V_inv_X(sigma, Eigen::VectorXd b, int n, int max_iteration, double tol);
        int n = eta.size();
        Eigen::VectorXd sigma_inv_eta = compute_V_inv_X(sigma, eta, n, max_iteration, tol);
        //Compute sigma_inv_X
        Eigen::MatrixXd sigma_inv_X = Eigen::MatrixXd::Zero(X.rows(),X.cols());
        for(int X_col=0;X_col<X.cols();X_col++){
            Eigen::VectorXd b_i = X.col(X_col);
            Eigen::VectorXd x = compute_V_inv_X(sigma, b_i, X.rows(), max_iteration, tol);
            sigma_inv_X.col(X_col) = x;
        }
        Eigen::SparseMatrix<double> X_T_sigma_inv_X = (X.transpose() * sigma_inv_X).sparseView();
        Eigen::VectorXd X_T_sigma_inv_eta = X.transpose() * sigma_inv_eta;
        beta_hat = compute_V_inv_X(X_T_sigma_inv_X, X_T_sigma_inv_eta, n, max_iteration, tol) ;
        std::cout << "beta_hat: " << beta_hat.transpose() << std::endl;

        // Update u_hat
        Eigen::VectorXd residual = eta - X * beta_hat;
        Eigen::VectorXd sigma_inv_residual = compute_V_inv_X(sigma, (eta - X*beta_hat), n, max_iteration, tol);
        u_hat = tau * K * sigma_inv_residual;
        cout<<"u_hat: "<<u_hat.head(10).transpose()<<endl;
        // Convergence criteria
        cout<<pre_tau<<"\t"<<tau<<endl;

        if ((std::abs(pre_tau - tau) < 1e-3)||(tau_next>10)) {
            std::cout << "Converged\n";
            std::cout << "Final tau: " << tau << std::endl;
            tau_ = tau;
            pi = sigmoid(X * beta_hat + u_hat);
            mu = pi.array() * total.array();
            eta = (X * beta_hat).array() + u_hat.array() + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));
            u_hat_ = u_hat;
            cout<<"pi: "<<pi.head(10).transpose()<<endl;
            cout<<"mu: "<<mu.head(10).transpose()<<endl;
            cout<<"eta: "<<eta.head(10).transpose()<<endl;
            break;
             // Return the result
        }
        pre_tau = tau;
        pre_phi = phi;
    }

    //std::cout << "Not converged!" << std::endl;
    }

    // facilitate computation
    Eigen::VectorXd compute_V_inv_X(Eigen::SparseMatrix<double> result, Eigen::VectorXd b, int n, int max_iteration, double tol);
    
    double logAbsDeterminant(const Eigen::SparseMatrix<double>& V);

    void read_in_genotype(const std::string vcf_file, const std::vector<std::pair<std::string, int>> &positions, 
                      int windowsize, const std::vector<std::string> common_sample) {
    std::cout << "Start reading in Genotype!" << std::endl;
    //vector<vector<int>> g;
    //vector<string> chr_pos;
    // Open VCF file
    htsFile *vcf = bcf_open(vcf_file.c_str(), "r");
    if (!vcf) {
        std::cerr << "Error opening VCF file." << std::endl;
        return;
    }

    // Initialize VCF reader
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    if (!hdr) {
        std::cerr << "Error reading VCF header." << std::endl;
        bcf_close(vcf);
        return;
    }

    // Load the VCF index (you need to have the index file available, e.g., .csi or .tbi)
    hts_idx_t *idx = bcf_index_load(vcf_file.c_str());
    if (!idx) {
        std::cerr << "Error loading VCF index." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return;
    }

    // Prepare to read records
    bcf1_t *record = bcf_init();
    if (!record) {
        std::cerr << "Error initializing VCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return;
    }

    // Create sample index map
    std::unordered_map<std::string, int> sample_idx_map;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        sample_idx_map[hdr->samples[i]] = i;
    }

    // Iterate over the positions
    for (const auto &pos_info : positions) {
        string chrom = pos_info.first;
        int pos = pos_info.second;

        // Fetch records within the windowsize
        int start = (pos - windowsize > 0) ? (pos - windowsize) : 0;
        int end = pos + windowsize;
        stringstream ss;
        ss<<chrom<<":"<<start<<"-"<<end;
        cout<<ss.str()<<endl;
        hts_itr_t *itr = bcf_itr_querys(idx, hdr, ss.str().c_str());
        if (!itr) {
            std::cerr << "Error creating iterator for " << chrom << ":" << start << "-" << end << std::endl;
            cout<<"Error"<<endl;
            return;
        }

        cout<<"iterator created!"<<endl;
        if(bcf_itr_next(vcf,itr,record)<0){
            cout<<"No records here"<<endl;
            bcf_destroy(record);
            bcf_hdr_destroy(hdr);
            bcf_close(vcf);
            return;
        }
        while (bcf_itr_next(vcf, itr, record) >= 0) {
            bcf_unpack(record, BCF_UN_ALL);  // Unpack record
            int pos = record->pos + 1;
            string ref=record->d.allele[0];
            string alt=record->d.allele[1];
            string snp = chrom + ":" + to_string(pos)+":"+ref+":"+alt;
            chr_pos.push_back(snp);
            // Add record ID to chr_pos
            Eigen::VectorXd tmp_g = Eigen::VectorXd::Zero(common_sample.size());

            // Iterate through common_sample and extract genotype dosage
            for (int c=0;c<common_sample.size();c++) {
                string sample = common_sample[c];
                std::string sample_substr = sample.substr(0, sample.find(":")); // Extract sample substring
                if (sample_idx_map.find(sample_substr) != sample_idx_map.end()) {
                    int sample_index = sample_idx_map[sample_substr];
                    int *gt_arr = nullptr, n_gts = 0;

                    // Extract genotype information
                    if (bcf_get_genotypes(hdr, record, &gt_arr, &n_gts) >= 0) {
                        int dosage = 0;

                        // Sum alleles to get dosage (ignoring missing values -1)
                        for (int j = sample_index * 2; j < sample_index * 2 + 2 && j < n_gts; ++j) {
                            if (gt_arr[j] != bcf_gt_missing) {
                                dosage += bcf_gt_allele(gt_arr[j]);
                            }
                        }
                        tmp_g[c] = dosage;
                    }

                    if (gt_arr) {
                        free(gt_arr);  // Free memory for genotype array
                    }
                }
            }

            // Append the dosage vector to g
            g.push_back(tmp_g);
        }

        hts_itr_destroy(itr);  // Destroy iterator
    }
    cout<<g[0].size()<<"\t"<<g.size()<<"\t"<<g[0][1]<<"\t"<<g[1][1]<<"\t"<<chr_pos.size()<<endl;
    cout<<"Genotype read in"<<endl;
    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}


void read_in_genotype_trans(const std::string vcf_file, const std::vector<string> &variant_id, 
                      const std::vector<std::string> common_sample,const std::string chrom) {
    std::cout << "Start reading in Genotype!" << std::endl;
    //vector<vector<int>> g;
    //vector<string> chr_pos;
    // Open VCF file
    htsFile *vcf = bcf_open(vcf_file.c_str(), "r");
    if (!vcf) {
        std::cerr << "Error opening VCF file." << std::endl;
        return;
    }

    // Initialize VCF reader
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    if (!hdr) {
        std::cerr << "Error reading VCF header." << std::endl;
        bcf_close(vcf);
        return;
    }

    // Load the VCF index (you need to have the index file available, e.g., .csi or .tbi)
    hts_idx_t *idx = bcf_index_load(vcf_file.c_str());
    if (!idx) {
        std::cerr << "Error loading VCF index." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return;
    }

    // Prepare to read records
    bcf1_t *record = bcf_init();
    if (!record) {
        std::cerr << "Error initializing VCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return;
    }

    // Create sample index map
    std::unordered_map<std::string, int> sample_idx_map;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        sample_idx_map[hdr->samples[i]] = i;
    }

    // Iterate over the positions
        // Fetch records within the windowsize
    while (bcf_read(vcf, hdr, record) == 0) {
            bcf_unpack(record, BCF_UN_ALL);  // Unpack record
            if (record->d.id && std::find(variant_id.begin(), variant_id.end(), record->d.id) != variant_id.end()){
            int pos = record->pos + 1;
            string ref=record->d.allele[0];
            string alt=record->d.allele[1];
            string snp = chrom + ":" + to_string(pos)+":"+ref+":"+alt;
            chr_pos.push_back(snp);
            // Add record ID to chr_pos
            Eigen::VectorXd tmp_g = Eigen::VectorXd::Zero(common_sample.size());

            // Iterate through common_sample and extract genotype dosage
            for (int c=0;c<common_sample.size();c++) {
                string sample = common_sample[c];
                std::string sample_substr = sample.substr(0, sample.find(":")); // Extract sample substring
                if (sample_idx_map.find(sample_substr) != sample_idx_map.end()) {
                    int sample_index = sample_idx_map[sample_substr];
                    int *gt_arr = nullptr, n_gts = 0;

                    // Extract genotype information
                    if (bcf_get_genotypes(hdr, record, &gt_arr, &n_gts) >= 0) {
                        int dosage = 0;

                        // Sum alleles to get dosage (ignoring missing values -1)
                        for (int j = sample_index * 2; j < sample_index * 2 + 2 && j < n_gts; ++j) {
                            if (gt_arr[j] != bcf_gt_missing) {
                                dosage += bcf_gt_allele(gt_arr[j]);
                            }
                        }
                        tmp_g[c] = dosage;
                    }

                    if (gt_arr) {
                        free(gt_arr);  // Free memory for genotype array
                    }
                }
            }

            // Append the dosage vector to g
            g.push_back(tmp_g);
        }

        // Destroy iterator
    }
    cout<<g[0].size()<<"\t"<<g.size()<<"\t"<<g[0][1]<<"\t"<<g[1][1]<<"\t"<<chr_pos.size()<<endl;
    cout<<"Genotype read in"<<endl;
    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
}

    double calculate_stddev(const std::vector<double>& data) {
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
        double variance = 0.0;
        for (double val : data) {
            variance += (val - mean) * (val - mean);
            }
        variance /= (data.size() - 1);  // Using Bessel's correction
        return std::sqrt(variance);
    }


    double dispersion_estimate(int times, Eigen::MatrixXd covariate_adjusted_geno, Eigen::VectorXd residuals, Eigen::VectorXd pi,Eigen::SparseMatrix<double> Identity){
        vector<double> test_statistics;
        vector<double> pvalue;
        vector<double> vari_total;
        const double maf = 0.5;
        random_device rd;
        mt19937 gen(rd());
        //uniform_int_distribution<> dis(0, 2); // Uniform distribution between 0 and 2
        binomial_distribution<> geno_dist(2,maf);
        double score_vector, info_matrix, test_statistic, variance;
        Eigen::VectorXd genotype, g;
        Eigen::VectorXd orig = Eigen::VectorXd::NullaryExpr(Identity.cols(), [&]() { return geno_dist(gen); });
        double observed=0;
        double expected=0;
        for(int i=0;i<times;i++){
            test_statistics.clear();
            for(int j=0;j<100;j++){
                // Simulate random genotypes
            //g = Eigen::VectorXd::NullaryExpr(Identity.cols(), [&]() { return geno_dist(gen); });
            static thread_local std::mt19937 gen(std::random_device{}());
            std::shuffle(orig.data(), orig.data() + orig.size(), gen);
            g=orig;
            genotype=Identity*g;
            // Adjust genotype
            g = genotype - covariate_adjusted_geno * genotype;
            //g = genotype;
            genotype = g.array() - g.mean();
            
            // Compute score vector and info matrix
            score_vector = (residuals.array() * genotype.array()).sum();
            info_matrix = (genotype.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();
            
            // Test statistic
            test_statistic = score_vector / std::sqrt(info_matrix);
            test_statistics.push_back(test_statistic);
            observed=observed+score_vector*score_vector;
	    expected=expected+info_matrix;
	    }
            variance = calculate_stddev(test_statistics);
            vari_total.push_back(variance);
        }
        double vari_mean = std::accumulate(vari_total.begin(), vari_total.end(), 0.0) / vari_total.size();
        cout<<"normalize:"<<vari_mean<<endl;
	cout<<"normalize2:"<<std::sqrt(observed/expected)<<endl;
        return std::sqrt(observed/expected);
    }

    void compute(string output_path, string site,Eigen::SparseMatrix<double> Identity,int times){  
        //Write in middle terms for pvalue computation: residuals, pi, total, y, dispersion
        //need eta, X, g, y, total
        //Compute covariate adjusted geno
        string file = output_path + "/" + site + ".middle";
        ofstream fout(file);
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<"\t"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
        Eigen::VectorXd pi = sigmoid(X*beta_hat+u_hat_);
        Eigen::VectorXd mu = pi.array()*total_.array();
        Eigen::VectorXd residuals = y_.array() - mu.array();
        cout<<"Start normalize parameter estimation"<<endl;
        double dispersion = dispersion_estimate(times,covariate_adjusted_geno,residuals,pi,Identity);
        fout<<"Splice_site"<<"\t"<<site<<"\t"<<dispersion<<"\t"<<tau_<<endl;
        fout<<"residuals"<<"\t";
        for(int i = 0;i<residuals.size();i++){
            fout<<residuals[i]<<"\t";
        }
        fout<<endl;
        fout<<"pi"<<"\t";
        for(int i = 0;i<pi.size();i++){
            fout<<pi[i]<<"\t";
        }
        fout<<endl;
        fout<<"total"<<"\t";
        for(int i = 0;i<total_.size();i++){
            fout<<total_[i]<<"\t";
        }
        fout<<endl;
        fout<<"y"<<"\t";
        for(int i = 0;i<y_.size();i++){
            fout<<y_[i]<<"\t";
        }
        fout<<endl;
        fout.close();
    }

   double computeRegressionSlope(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    // Ensure the vectors are of the same size
    if (x.size() != y.size()) {
        throw std::invalid_argument("Vectors x and y must have the same size.");
    }

    // Calculate the mean of x and y
    double x_mean = x.mean();
    double y_mean = y.mean();

    // Calculate the standard deviation of x and y
    double x_stddev = std::sqrt((x.array() - x_mean).square().mean());
    double y_stddev = std::sqrt((y.array() - y_mean).square().mean());

    // Normalize x and y
    Eigen::VectorXd x_normalized = (x.array() - x_mean) / x_stddev;
    Eigen::VectorXd y_normalized = (y.array() - y_mean) / y_stddev;

    // Compute the numerator and denominator for the slope
    double numerator = x_normalized.dot(y_normalized);
    double denominator = x_normalized.squaredNorm();

    // Calculate and return the slope
    return numerator / denominator;
    }

    void pvalue_beta_sd_compute(string output_file, string site, double dispersion, double threshold, Eigen::VectorXd residuals,Eigen::VectorXd pi,Eigen::VectorXd y,Eigen::VectorXd total){
        if(chr_pos.size()==0){
            cout<<"No genotype here"<<endl;
            return;
        }
        cout<<"Start pvalue computing!"<<endl;
        ofstream fout(output_file);
        Eigen::VectorXd tmp_t_p_1_p = total_.array() * pi.array() * (1.0 - pi.array());
        Eigen::VectorXd tmp_p_1_p = pi.array() * (1.0 - pi.array());
        int sample_size=residuals.size();
        Eigen::VectorXd pheno = residuals.array()/total_.array();
        //covariate_adjusted_genotype
        Eigen::VectorXd W_vec = ((pi.array()*(1-pi.array()))).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<"\t"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
        //
        for(int i=0;i<chr_pos.size();i++){
            const auto &g_vec = g[i];
    // Copy values from std::vector to Eigen::VectorXd
            Eigen::VectorXd geno = g_vec - covariate_adjusted_geno * g_vec;
            //Eigen::VectorXd geno = g_vec;
            Eigen::VectorXd genotype = geno.array() - geno.mean();
            
            #pragma omp parallel
            {
            // Compute score vector and info matrix
            double score_vector = (residuals.array() * genotype.array()).sum();
            //double info_matrix = (genotype.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();
            double info_matrix = (genotype.array().square() * tmp_t_p_1_p.array() ).sum();

            // Test statistic
            double test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
            double p_value = gsl_cdf_chisq_Q(test_statistics*test_statistics,1);

    // 2. Compute variance of test statistics
            double variance_test = std::sqrt((genotype.array().square() * tmp_t_p_1_p.array() ).sum()) * dispersion;

    // 3. Compute effect size
            //double effect = test_statistics / variance_test;
            //double effect = test_statistics/std::sqrt(sample_size);
    // 4. Compute z_score using the inverse of the normal CDF (ppf in Python)
            double z_score = gsl_cdf_gaussian_Pinv(1.0 - p_value / 2, 1.0);
            double effect = computeRegressionSlope(genotype,pheno);
    // 5. Compute standard error
            //double standard_error = std::abs(effect / test_statistics);
            double standard_error = std::abs(effect/z_score);
            #pragma imp critical
            {
            if(p_value<threshold){
                fout<<site<<"\t"<<chr_pos[i]<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;}
            }
            //fout<<site<<"\t"<<chr_pos[i]<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;}
            }
        }
            fout.close();
    }

    void DS(string output_file, string site, double dispersion, vector<int> group, Eigen::VectorXd residuals,Eigen::VectorXd pi,Eigen::VectorXd y,Eigen::VectorXd total){
        ofstream fout(output_file);
        Eigen::VectorXd tmp_t_p_1_p = total_.array() * pi.array() * (1.0 - pi.array());
        Eigen::VectorXd tmp_p_1_p = pi.array() * (1.0 - pi.array());
        Eigen::VectorXd group_double(group.size());
        Eigen::VectorXd pheno = residuals.array()/total_.array();
        for(size_t i=0;i<group.size();i++){
            group_double[i]=static_cast<double>(group[i]);
        }
        Eigen::VectorXd W_vec = ((pi.array()*(1-pi.array()))).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<"\t"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
    // Copy values from std::vector to Eigen::VectorXd
        Eigen::VectorXd geno_double = group_double - covariate_adjusted_geno * group_double;
        Eigen::VectorXd group_new = geno_double.array() - geno_double.mean();
            // Compute score vector and info matrix
        double score_vector = (residuals.array() * group_new.array()).sum();
            //double info_matrix = (genotype.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();
        double info_matrix = (group_new.array().square() * tmp_t_p_1_p.array() ).sum();

            // Test statistic
        double test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
        double p_value = gsl_cdf_chisq_Q(test_statistics*test_statistics,1);

    // 2. Compute variance of test statistics
        double variance_test = std::sqrt((group_new.array().square() * tmp_p_1_p.array() ).sum()) * dispersion;

    // 3. Compute effect size
        double effect = computeRegressionSlope(group_new,pheno);

    // 4. Compute z_score using the inverse of the normal CDF (ppf in Python)
        double z_score = gsl_cdf_gaussian_Pinv(1.0 - p_value / 2, 1.0);

    // 5. Compute standard error
        double standard_error = std::abs(effect / z_score);

        fout<<site<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;
        fout.close();
    }


string get_info_as_string(bcf_hdr_t* hdr, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_INFO);

    std::ostringstream infoString;

    for (int i = 0; i < record->n_info; ++i) {
        bcf_info_t* info = &record->d.info[i];
        const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);

        if (i > 0) infoString << ";"; // separate multiple INFO fields

        infoString << key << "=";

        if (info->len == 0) continue;

        // Handle different data types
        if (info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32) {
            int32_t* values = (int32_t*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                infoString << values[j];
            }
        } else if (info->type == BCF_BT_FLOAT) {
            float* values = (float*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                infoString << values[j];
            }
        } else if (info->type == BCF_BT_CHAR) {
            infoString << (char*)info->vptr;
        }
    }

    return infoString.str();
}


int read_in_genotype_rare(const std::string vcf_file, const std::vector<string> &annotation, const std::vector<std::pair<std::string, int>> &positions,
                      const std::vector<std::string> common_sample,const std::string chrom, int windowsize) {
    std::cout << "Start reading in Genotype!" << std::endl;
    //vector<vector<int>> g;
    //vector<string> chr_pos;
    // Open VCF file
    htsFile *vcf = bcf_open(vcf_file.c_str(), "r");
    if (!vcf) {
        std::cerr << "Error opening VCF file." << std::endl;
        return -1;
    }

    // Initialize VCF reader
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    if (!hdr) {
        std::cerr << "Error reading VCF header." << std::endl;
        bcf_close(vcf);
        return -1;
    }

    // Load the VCF index (you need to have the index file available, e.g., .csi or .tbi)
    hts_idx_t *idx = bcf_index_load(vcf_file.c_str());
    if (!idx) {
        std::cerr << "Error loading VCF index." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return -1;
    }

    // Prepare to read records
    bcf1_t *record = bcf_init();
    if (!record) {
        std::cerr << "Error initializing VCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return -1;
    }

   // Create sample index map
    std::unordered_map<std::string, int> sample_idx_map;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        sample_idx_map[hdr->samples[i]] = i;
    }

    // Iterate over the positions
    std::pair<std::string, int> pos_info = positions[0];
    string chr = pos_info.first;
    int pos = pos_info.second;

        // Fetch records within the windowsize
    int start = (pos - windowsize > 0) ? (pos - windowsize) : 0;
    int end = pos + windowsize;
    stringstream ss;
    ss<<chrom<<":"<<start<<"-"<<end;
    cout<<ss.str()<<endl;
    hts_itr_t *itr = bcf_itr_querys(idx, hdr, ss.str().c_str());
    if (!itr) {
            std::cerr << "Error creating iterator for " << chrom << ":" << start << "-" << end << std::endl;
            cout<<"Error"<<endl;
            return -1;
    }

    cout<<"iterator created!"<<endl;
    if(bcf_itr_next(vcf,itr,record)<0){
            cout<<"No records here"<<endl;
            bcf_destroy(record);
            bcf_hdr_destroy(hdr);
            bcf_close(vcf);
            return -1;
    }
    //define two vectors, MAF & distance
    double maf,dist;
    int pos_site=pos;
    while (bcf_itr_next(vcf, itr, record) >= 0) {
            bcf_unpack(record, BCF_UN_ALL);  // Unpack record
            int pos = record->pos + 1;
            dist = abs(pos-pos_site);
            distance.push_back(dist);
            string ref=record->d.allele[0];
            string alt=record->d.allele[1];
            // obtain variant annotation information 
            string variant_annot = get_info_as_string(hdr, record);
            string maf_tmp=variant_annot.substr(variant_annot.find("MAF"));
            string maf_tmp_2=maf_tmp.substr(0,maf_tmp.find(";"));
            maf=stod(maf_tmp_2.substr(maf_tmp_2.find("=")+1));
            MAF.push_back(maf);

            //bool exist_variant=false;
            //for(int e=0;e<annotation.size();e++){
            //    if(variant_annot.find(annotation[e])!=string::npos){
            //        exist_variant=true;
            //    }
            //}

            //need to maintain a vector with MAF and distance to splice site here 2025.12.8
            bool exist_variant=true;
            if(exist_variant==true){
            string snp = chrom + ":" + to_string(pos)+":"+ref+":"+alt;
            chr_pos.push_back(snp);
            cout<<snp<<endl;
            // Add record ID to chr_pos
            Eigen::VectorXd tmp_g = Eigen::VectorXd::Zero(common_sample.size());

            // Iterate through common_sample and extract genotype dosage
            for (int c=0;c<common_sample.size();c++) {
                string sample = common_sample[c];
                std::string sample_substr = sample.substr(0, sample.find(":")); // Extract sample substring
                if (sample_idx_map.find(sample_substr) != sample_idx_map.end()) {
                    int sample_index = sample_idx_map[sample_substr];
                    int *gt_arr = nullptr, n_gts = 0;

                    // Extract genotype information
                    if (bcf_get_genotypes(hdr, record, &gt_arr, &n_gts) >= 0) {
                        int dosage = 0;

                        // Sum alleles to get dosage (ignoring missing values -1)
                        for (int j = sample_index * 2; j < sample_index * 2 + 2 && j < n_gts; ++j) {
                            if (gt_arr[j] != bcf_gt_missing) {
                                dosage += bcf_gt_allele(gt_arr[j]);
                            }
                        }
                        tmp_g[c] = dosage;
                    }

                    if (gt_arr) {
                        free(gt_arr);  // Free memory for genotype array
                    }
                }
            }

            // Append the dosage vector to g
            g.push_back(tmp_g);
    }}
    hts_itr_destroy(itr);  // Destroy iterator
    cout<<"destroy itr"<<endl;
    
    cout<<g[0].size()<<"\t"<<g.size()<<"\t"<<g[0][1]<<"\t"<<chr_pos.size()<<endl;
    cout<<"Genotype read in"<<endl;
    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);
    return g.size();
}

double beta_pdf(double x, double alpha, double beta) {
    if (x < 0.0 || x > 1.0) {
        return 0.0; // Beta distribution is defined on the interval [0, 1]
    }

    // Compute Beta function B(alpha, beta) = Gamma(alpha) * Gamma(beta) / Gamma(alpha + beta)
    double beta_function = std::tgamma(alpha) * std::tgamma(beta) / std::tgamma(alpha + beta);

    // Beta PDF formula
    return (std::pow(x, alpha - 1) * std::pow(1 - x, beta - 1)) / beta_function;
}

void collapse_genotype(int windowsize){
    // incorporate both MAF & weight
    vector<Eigen::VectorXd> g_tmp = g;
    int n = g.front().size();
    Eigen::VectorXd result = Eigen::VectorXd::Constant(n,0);
    for(int i=0;i<distance.size();i++){
        //double weight_MAF = beta_pdf(MAF[i], alpha, beta);;
        //double weight_distance = beta_pdf(distance[i]/windowsize, alpha, beta);;
        double weight= exp(-(distance[i]*MAF[i]));
        //cout<<"weight"<<weight_MAF<<"\t"<<weight_distance<<endl;
        for(int j=0;j<g[i].size();j++){
            g[i][j]*=weight;
        }
        result=result+g[i];
    }
    //for (const auto& v : g) {
    //    result = result.cwiseMax(v); // elementwise max
    //}
    g.clear();
    /*double minVal=result.minCoeff();
    double maxVal=result.maxCoeff();
    Eigen::VectorXd scaled = (result.array() - minVal) / (maxVal - minVal);
    for (int i = 0; i < scaled.size(); ++i) {
        if (scaled[i] < 0.33) {
            scaled[i] = 0;
        } else if (scaled[i] < 0.66) {
            scaled[i] = 1;
        } else {
            scaled[i] = 2;
        }
    }*/
    //g.push_back(scaled);
    g.push_back(result);
    Eigen::VectorXd result_tmp=result;
    double weight;
    for(int i=0;i<distance.size();i++){
        //double weight_MAF = beta_pdf(MAF[i], alpha, beta);;
        //double weight_distance = beta_pdf(distance[i]/windowsize, alpha, beta);;
        weight= exp(-(distance[i]*MAF[i]));
        //cout<<"weight"<<weight_MAF<<"\t"<<weight_distance<<endl;
        for(int j=0;j<g_tmp[i].size();j++){
            g_tmp[i][j]*=weight;
        }
        result_tmp=result-g_tmp[i];
        /*minVal=result_tmp.minCoeff();
        maxVal=result_tmp.maxCoeff();
        scaled = (result_tmp.array() - minVal) / (maxVal - minVal);
        for (int i = 0; i < scaled.size(); ++i) {
         if (scaled[i] < 0.33) {
            scaled[i] = 0;
         } else if (scaled[i] < 0.66) {
            scaled[i] = 1;
         } else {
            scaled[i] = 2;
        }
        }*/
        g.push_back(result_tmp);
    }
    cout<<"genotype collapse"<<endl;
    cout<<g[0].size()<<"\t"<<g.size()<<"\t"<<g[0][1]<<"\t"<<chr_pos.size()<<endl;
}

double compute_MAF(const Eigen::VectorXd& geno) {
    // Count the occurrences of 0s, 1s, and 2s in geno
    double count_0 = (geno.array() == 0).count(); // number of homozygous reference
    //double count_1 = (geno.array() == 1).count(); // number of heterozygous
    //double count_2 = (geno.array() == 2).count(); // number of homozygous alternate
    double count_2 = (geno.array() != 0).count();
    // Total number of alleles (2 alleles per genotype)
    double total_alleles = 2 * geno.size();
    
    // Count of reference (0) and alternate (2) alleles
    double count_reference_allele = count_0 * 2 ; // Each homozygous reference contributes 2, heterozygous contributes 1
    double count_alternate_allele = count_2 * 2 ; // Each homozygous alternate contributes 2, heterozygous contributes 1

    // MAF is the frequency of the minor allele (the less frequent one)
    double maf = std::min(count_reference_allele, count_alternate_allele) / total_alleles;
    
    return maf;
}

void pvalue_beta_sd_compute_rare(string output_file, string site, double dispersion, double threshold, Eigen::VectorXd residuals,Eigen::VectorXd pi,Eigen::VectorXd y,Eigen::VectorXd total){
        if(chr_pos.size()==0){
            cout<<"No genotype here"<<endl;
            return;
        }
        cout<<"Start pvalue computing!"<<endl;
        ofstream fout(output_file);
        Eigen::VectorXd tmp_t_p_1_p = total_.array() * pi.array() * (1.0 - pi.array());
        Eigen::VectorXd tmp_p_1_p = pi.array() * (1.0 - pi.array());
        int sample_size=residuals.size();
        Eigen::VectorXd pheno = residuals.array()/total_.array();
        //covariate_adjusted_genotype
        Eigen::VectorXd W_vec = ((pi.array()*(1-pi.array()))).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<"\t"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
        //
        double maf = 0;
        const auto &g_vec = g[0];
        maf = compute_MAF(g_vec);
    // Copy values from std::vector to Eigen::VectorXd
        Eigen::VectorXd geno = g_vec - covariate_adjusted_geno * g_vec;
            //Eigen::VectorXd geno = g_vec;
        Eigen::VectorXd genotype = geno.array() - geno.mean();
            // Compute score vector and info matrix
        double score_vector = (residuals.array() * genotype.array()).sum();
            //double info_matrix = (genotype.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();
        double info_matrix = (genotype.array().square() * tmp_t_p_1_p.array() ).sum();

            // Test statistic
        double test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
        double p_value = gsl_cdf_chisq_Q(test_statistics*test_statistics,1);

    // 2. Compute variance of test statistics
        double variance_test = std::sqrt((genotype.array().square() * tmp_t_p_1_p.array() ).sum()) * dispersion;

    // 3. Compute effect size
            //double effect = test_statistics / variance_test;
            //double effect = test_statistics/std::sqrt(sample_size);
    // 4. Compute z_score using the inverse of the normal CDF (ppf in Python)
        double z_score = gsl_cdf_gaussian_Pinv(1.0 - p_value / 2, 1.0);
        double effect = computeRegressionSlope(genotype,pheno);
    // 5. Compute standard error
            //double standard_error = std::abs(effect / test_statistics);
        double standard_error = std::abs(effect/z_score);
        #pragma imp critical
        {
            if(p_value<threshold){
                fout<<site<<"\t"<<"burden"<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<"\t"<<maf<<endl;}
        }
            //fout<<site<<"\t"<<chr_pos[i]<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;}

        double total_test_statistics=test_statistics;
        double relative_contribution;
        for(int i=0;i<chr_pos.size();i++){
            const auto &g_vec = g[i+1];
            maf = compute_MAF(g_vec);
    // Copy values from std::vector to Eigen::VectorXd
            geno = g_vec - covariate_adjusted_geno * g_vec;
            //Eigen::VectorXd geno = g_vec;
            genotype = geno.array() - geno.mean();
            
            #pragma omp parallel
            {
            // Compute score vector and info matrix
            score_vector = (residuals.array() * genotype.array()).sum();
            //double info_matrix = (genotype.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();
            info_matrix = (genotype.array().square() * tmp_t_p_1_p.array()).sum();

            // Test statistic
            test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
            relative_contribution = 1-(test_statistics/total_test_statistics)*(test_statistics/total_test_statistics);
    // 5. Compute standard error
            //double standard_error = std::abs(effect / test_statistics);
            #pragma imp critical
            {
            if(p_value<threshold){
                fout<<site<<"\t"<<chr_pos[i]<<"\t"<<relative_contribution<<"\t"<<maf<<endl;}
            }
            //fout<<site<<"\t"<<chr_pos[i]<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;}
            }
        }
            fout.close();
    }

private:
    Eigen::VectorXd& eta;
    Eigen::VectorXd& u_hat_;
    const Eigen::MatrixXd& X;
    Eigen::VectorXd& beta_hat;
    const Eigen::SparseMatrix<double> K;
    Eigen::SparseMatrix<double> W_;
    Eigen::MatrixXd& sigma_inv_X;
    const int max_iteration;
    const double tol;
    vector<Eigen::VectorXd> g;
    vector<string> chr_pos;
    vector<double> MAF;
    vector<double> distance;
    Eigen::VectorXd& y_;
    Eigen::VectorXd& total_;
    double tau_;
};

int model(int argc, char *argv[]);
int QTL_parse_options(int argc, char *argv[]);
int DS_parse_options(int argc, char *argv[]);
int QTL(int argc, char *argv[]);
int trans_QTL(int argc, char *argv[]);
int DS(int argc, char *argv[]);
int rare_mapping(int argc, char *argv[]);
int rare_mapping_parse_options(int argc, char *argv[]);


#endif

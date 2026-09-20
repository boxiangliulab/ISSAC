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
        double GRM_val_;
        // iteration times for estimation fixed effect and random effect
        int iter_time_;
        double iter_thre_;
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
                  Eigen::VectorXd& total_, Eigen::SparseMatrix<double> K, const int max_iteration_,
                  const double tol_, Eigen::MatrixXd& sigma_inv_X, double tau_) 
        : eta(eta), X(X), beta_hat(beta_hat), u_hat_(u_hat_),y_(y_),total_(total_), K(K), max_iteration_(max_iteration_), tol_(tol_), sigma_inv_X(sigma_inv_X), tau_(tau_) {}

static double reml_criterion(const vector<double> &theta, vector<double> &grad, void *data){
    auto *args = static_cast<std::tuple<Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd,
                                         Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>,
                                         Eigen::SparseMatrix<double>> *>(data);
    Eigen::VectorXd eta            = get<0>(*args);
    Eigen::MatrixXd X              = get<1>(*args);
    Eigen::VectorXd beta_hat       = get<2>(*args);
    Eigen::SparseMatrix<double> W_inv = get<3>(*args);
    Eigen::SparseMatrix<double> K     = get<4>(*args);
    Eigen::SparseMatrix<double> I     = get<5>(*args);

    double tau_g = std::exp(theta[0]);
    double tau_o = std::exp(theta[1]);

    Eigen::SparseMatrix<double> result = W_inv + tau_g * K + tau_o * I;

    // Factorize V once, reuse for all solves below
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(result);
    solver.factorize(result);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition of V failed at tau_g=" << tau_g
                   << ", tau_o=" << tau_o << std::endl;
        return std::numeric_limits<double>::max();
    }
    double log_det_V = solver.logAbsDeterminant();

    // Solve V^{-1} X column by column, reusing the same factorization
    Eigen::MatrixXd V_inv_X = Eigen::MatrixXd::Zero(X.rows(), X.cols());
    for (int i = 0; i < X.cols(); i++) {
        V_inv_X.col(i) = solver.solve(X.col(i));
    }

    Eigen::SparseMatrix<double> X_t_V_inv_X = (X.transpose() * V_inv_X).sparseView();

    Eigen::VectorXd residual = eta - X * beta_hat;
    Eigen::VectorXd V_inv_residual = solver.solve(residual);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_2;
    solver_2.analyzePattern(X_t_V_inv_X);
    solver_2.factorize(X_t_V_inv_X);
    if (solver_2.info() != Eigen::Success) {
        std::cerr << "Decomposition of X'V^{-1}X failed at tau_g=" << tau_g
                   << ", tau_o=" << tau_o << std::endl;
        return std::numeric_limits<double>::max();
    }
    double log_det_X_t_V_inv_X = solver_2.logAbsDeterminant();

    double reml_value = -0.5 * (log_det_V + log_det_X_t_V_inv_X
                                 + residual.transpose() * V_inv_residual);
    return -reml_value;
}

    Eigen::VectorXd sigmoid(const Eigen::VectorXd & x){
        return 1.0/(1.0+(-x.array()).exp());
    }

void update(Eigen::VectorXd& y, Eigen::VectorXd& total){

    y_ = y;
    total_ = total;

    double tau_g = 0;
    double tau_o = 0;
    convergence_status_ = false;
    Eigen::VectorXd start_beta = beta_hat; 
    
    Eigen::VectorXd pre_beta = beta_hat;
    double pre_tau_g = 0;
    double pre_tau_o = 0;

    // Fixed inner CG solver settings, used consistently for every V^{-1} * (.) solve below
    const int    CG_MAX_ITERATION = 1000;
    const double CG_TOL = 1e-6;

    int n = y.size();
    Eigen::VectorXd u_hat_g = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd u_hat_o = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd u_hat   = Eigen::VectorXd::Zero(n); // u_hat_g + u_hat_o

    Eigen::VectorXd pi = sigmoid(X * beta_hat + u_hat);
    Eigen::VectorXd mu = pi.array() * total.array();

    Eigen::SparseMatrix<double> I_mat(n, n);
    I_mat.setIdentity();

    cout << "Start iteration" << endl;
    Eigen::VectorXd W_vec,W_inv_vec;
    Eigen::MatrixXd W_diag,W_inv_diag;
    Eigen::SparseMatrix<double> W_inv_sparse;
    for (int i = 0; i < max_iteration_; ++i) {
        pi = sigmoid(X * beta_hat + u_hat);
        mu = pi.array() * total.array();
        eta = (X * beta_hat).array() + u_hat.array()
              + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));

        W_vec = (total.array() * (pi.array()*(1-pi.array()))).matrix();
        W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();

        W_inv_vec = (1/(total.array() * (pi.array()*(1-pi.array())))).matrix();
        W_inv_diag = W_inv_vec.asDiagonal();
        W_inv_sparse = W_inv_diag.sparseView();

        // initialize theta = [log(tau_g), log(tau_o)]
        double mean_eta = eta.mean();
        double variance_eta = (eta.array() - mean_eta).square().mean();
        double log_var_eta = log(variance_eta);
        cout << "log_var_eta\t" << log_var_eta << endl;
        if (std::isnan(log_var_eta)) {
            throw std::invalid_argument("Iteration failed");
        }

        vector<double> theta(2, log_var_eta - log(2.0));

        cout << "eta" << "\t" << eta.head(10).transpose() << endl;
        cout << "pi" << "\t" << pi.head(10).transpose() << endl;
        cout << "y" << "\t" << y.head(10).transpose() << endl;
        cout << "total" << "\t" << total.head(10).transpose() << endl;
        cout << "Defined successful" << endl;

        auto args = std::make_tuple(eta, X, beta_hat, W_inv_sparse, K, I_mat);
        nlopt::opt opt(nlopt::LN_BOBYQA, 2);
        std::vector<double> lb(2, std::log(1e-6));
        std::vector<double> ub(2, std::log(100));
        opt.set_lower_bounds(lb);
        opt.set_upper_bounds(ub);
        opt.set_min_objective(reml_criterion, &args);
        opt.set_xtol_rel(1e-4);    

        double minf;
        bool opt_success = true;
        try {
            nlopt::result result = opt.optimize(theta, minf);
            std::cout << "Found minimum at f(x) = " << minf << std::endl;
            std::cout << "Optimal parameters: log_tau_g=" << theta[0] << std::endl;
        } catch (std::exception &e) {
            std::cerr << "NLopt failed at iteration " << i << ": " << e.what() << std::endl;
            opt_success = false;
        }

        if (!opt_success) {
            std::cerr << "Warning: variance component optimization failed at iteration " << i
                       << "; treating as non-convergence." << std::endl;
            convergence_status_ = false;
            break;
        }

        double tau_g_next = std::exp(theta[0]);
        double tau_o_next = std::exp(theta[1]);

        std::cout << "Iteration " << i << " : tau=" << tau_g_next << std::endl;
        
        const double divergence_threshold = 95.0;
        bool diverged = (tau_g_next >= divergence_threshold) || (tau_o_next >= divergence_threshold);
        if (diverged) {
        std::cerr << "Warning: variance component diverged at iteration " << i
                   << " (tau_g=" << tau_g_next << ")" << std::endl;
        convergence_status_ = false;
        break;}

        pre_tau_g = tau_g;
        pre_tau_o = tau_o;
        tau_g = tau_g_next;
        tau_o = tau_o_next;

        sigma  = W_inv_sparse + tau_g * K + tau_o * I_mat;

        pre_beta = beta_hat;
        int n_ = eta.size();
        Eigen::VectorXd sigma_inv_eta = compute_V_inv_X(sigma, eta, n_, CG_MAX_ITERATION, CG_TOL);

        sigma_inv_X = Eigen::MatrixXd::Zero(X.rows(), X.cols());
        for (int X_col = 0; X_col < X.cols(); X_col++) {
            Eigen::VectorXd b_i = X.col(X_col);
            Eigen::VectorXd x = compute_V_inv_X(sigma, b_i, X.rows(), CG_MAX_ITERATION, CG_TOL);
            sigma_inv_X.col(X_col) = x;
        }
        Eigen::SparseMatrix<double> X_T_sigma_inv_X = (X.transpose() * sigma_inv_X).sparseView();
        Eigen::VectorXd X_T_sigma_inv_eta = X.transpose() * sigma_inv_eta;
        beta_hat = compute_V_inv_X(X_T_sigma_inv_X, X_T_sigma_inv_eta, n_, CG_MAX_ITERATION, CG_TOL);
        std::cout << "beta_hat: " << beta_hat.transpose() << std::endl;

        Eigen::VectorXd residual = eta - X * beta_hat;
        Eigen::VectorXd sigma_inv_residual = compute_V_inv_X(sigma, residual, n_, CG_MAX_ITERATION, CG_TOL);
        u_hat_g = tau_g * (K * sigma_inv_residual);
        u_hat_o = tau_o * sigma_inv_residual; // I * sigma_inv_residual = sigma_inv_residual
        u_hat = u_hat_g + u_hat_o;

        cout << "u_hat: " << u_hat.head(10).transpose() << endl;
        //cout << "u_hat_o: " << u_hat_o.head(10).transpose() << endl;
        //cout << pre_tau_g << "\t" << tau_g << "\t" << pre_tau_o << "\t" << tau_o << endl;

        double beta_delta = (beta_hat - pre_beta).norm() / (pre_beta.norm() + 1e-10);
        bool tau_converged = (std::abs(pre_tau_g - tau_g) < tol_) && (std::abs(pre_tau_o - tau_o) < tol_);
        bool beta_converged = beta_delta < tol_;

        if (tau_converged && beta_converged) {
            std::cout << "Converged\n";
            std::cout << "Final tau_g: " << tau_g << std::endl;
            convergence_status_=true;
            tau_   = tau_g;   
            tau_o_ = tau_o;  
            pi = sigmoid(X * beta_hat + u_hat);
            mu = pi.array() * total.array();
            eta = (X * beta_hat).array() + u_hat.array()
                  + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));
            u_hat_ = u_hat;
            W_vec = (total.array() * (pi.array()*(1-pi.array()))).matrix();
            W_diag = W_vec.asDiagonal();
            W_ = W_diag.sparseView();
            cout << "pi: " << pi.head(10).transpose() << endl;
            cout << "mu: " << mu.head(10).transpose() << endl;
            cout << "eta: " << eta.head(10).transpose() << endl;
            break;
        }
    }

    //only binomial fixed model
    if(convergence_status_==false){
        u_hat   = Eigen::VectorXd::Zero(n);
        tau_ = 0.0;     
        tau_o_ = 0.0;
        beta_hat=start_beta;
        pi = sigmoid(X * beta_hat + u_hat);
        mu = pi.array() * total.array();
        eta = (X * beta_hat).array() + u_hat.array()
              + (y.array() - mu.array()) / (total.array() * pi.array() * (1.0 - pi.array()));
        u_hat_ = u_hat;
        W_vec = (total.array() * (pi.array()*(1-pi.array()))).matrix();
        W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        cout << "pi: " << pi.head(10).transpose() << endl;
        cout << "mu: " << mu.head(10).transpose() << endl;
        cout << "eta: " << eta.head(10).transpose() << endl;
    }
}

    // facilitate computation
    Eigen::VectorXd compute_V_inv_X(Eigen::SparseMatrix<double> result, Eigen::VectorXd b, int n, int max_iteration, double tol);
    
    double logAbsDeterminant(const Eigen::SparseMatrix<double>& V);

  bool read_in_genotype(
    const std::string vcf_file,
    const std::vector<std::pair<std::string, int>> &positions,
    int windowsize,
    const std::vector<std::string> common_sample
) {
    std::cout << "Start reading in Genotype!" << std::endl;

    // ============================================================
    // 1. Open VCF file
    // ============================================================
    htsFile *vcf = bcf_open(vcf_file.c_str(), "r");

    if (!vcf) {
        std::cerr << "Error: cannot open VCF file: "
                  << vcf_file << std::endl;
        return false;
    }


    // ============================================================
    // 2. Read VCF header
    // ============================================================
    bcf_hdr_t *hdr = bcf_hdr_read(vcf);

    if (!hdr) {
        std::cerr << "Error: cannot read VCF header." << std::endl;
        bcf_close(vcf);
        return false;
    }


    // ============================================================
    // 3. Load VCF index
    // ============================================================
    hts_idx_t *idx = bcf_index_load(vcf_file.c_str());

    if (!idx) {
        std::cerr << "Error: cannot load VCF index for: "
                  << vcf_file << std::endl;

        bcf_hdr_destroy(hdr);
        bcf_close(vcf);

        return false;
    }


    // ============================================================
    // 4. Initialize VCF record
    // ============================================================
    bcf1_t *record = bcf_init();

    if (!record) {
        std::cerr << "Error: cannot initialize VCF record."
                  << std::endl;

        hts_idx_destroy(idx);
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);

        return false;
    }


    // ============================================================
    // 5. Create sample index map
    // ============================================================
    std::unordered_map<std::string, int> sample_idx_map;

    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        sample_idx_map[hdr->samples[i]] = i;
    }


    // ============================================================
    // 6. Iterate over splice-site positions
    // ============================================================
    for (const auto &pos_info : positions) {

        std::string chrom = pos_info.first;
        int pos = pos_info.second;

        int start = (pos - windowsize > 0)
                        ? (pos - windowsize)
                        : 0;

        int end = pos + windowsize;

        std::stringstream ss;
        ss << chrom << ":" << start << "-" << end;

        std::cout << "Searching region: "
                  << ss.str() << std::endl;


        // --------------------------------------------------------
        // Create iterator for this region
        // --------------------------------------------------------
        hts_itr_t *itr =
            bcf_itr_querys(idx, hdr, ss.str().c_str());

        if (!itr) {
            std::cerr
                << "Warning: cannot create VCF iterator for "
                << chrom << ":"
                << start << "-"
                << end
                << ". Skipping this region."
                << std::endl;

            continue;
        }

        bool found_any = false;

        std::cout << "Iterator created!" << std::endl;


        // ========================================================
        // 7. Read variants in the window
        // ========================================================
        while (bcf_itr_next(vcf, itr, record) >= 0) {

            bcf_unpack(record, BCF_UN_ALL);

            // Need at least REF and ALT
            if (record->n_allele < 2) {
                std::cerr
                    << "Warning: record with fewer than 2 alleles "
                    << "encountered. Skipping."
                    << std::endl;
                continue;
            }

            int record_pos = record->pos + 1;

            std::string ref = record->d.allele[0];
            std::string alt = record->d.allele[1];

            std::string snp =
                chrom + ":" +
                std::to_string(record_pos) + ":" +
                ref + ":" +
                alt;


            // ----------------------------------------------------
            // Extract genotype data
            // ----------------------------------------------------
            int *gt_arr = nullptr;
            int n_gts = 0;

            int gt_call_result =
                bcf_get_genotypes(
                    hdr,
                    record,
                    &gt_arr,
                    &n_gts
                );


            // If genotype cannot be extracted, skip this variant
            if (gt_call_result < 0 || gt_arr == nullptr) {

                std::cerr
                    << "Warning: failed to extract genotypes for "
                    << snp
                    << ". Skipping this variant."
                    << std::endl;

                if (gt_arr) {
                    free(gt_arr);
                }

                continue;
            }


            // ----------------------------------------------------
            // Initialize genotype dosage vector
            // ----------------------------------------------------
            Eigen::VectorXd tmp_g =
                Eigen::VectorXd::Zero(common_sample.size());


            // ----------------------------------------------------
            // Extract genotype dosage for common samples
            // ----------------------------------------------------
            for (int c = 0;
                 c < static_cast<int>(common_sample.size());
                 ++c) {

                std::string sample = common_sample[c];

                std::string sample_substr =
                    sample.substr(0, sample.find(":"));

                auto sample_it =
                    sample_idx_map.find(sample_substr);

                if (sample_it == sample_idx_map.end()) {
                    continue;
                }

                int sample_index = sample_it->second;

                int dosage = 0;

                // Diploid genotype: two alleles per sample
                for (int j = sample_index * 2;
                     j < sample_index * 2 + 2 &&
                     j < n_gts;
                     ++j) {

                    if (gt_arr[j] != bcf_gt_missing) {
                        dosage +=
                            bcf_gt_allele(gt_arr[j]);
                    }
                }

                tmp_g[c] = dosage;
            }


            // ----------------------------------------------------
            // Free genotype array
            // ----------------------------------------------------
            free(gt_arr);


            // ----------------------------------------------------
            // Variant successfully read
            // ----------------------------------------------------
            chr_pos.push_back(snp);
            g.push_back(tmp_g);

            found_any = true;
        }


        // ========================================================
        // 8. No usable variants in this region
        // ========================================================
        if (!found_any) {

            std::cerr
                << "Warning: no usable genotype records found for "
                << chrom << ":"
                << start << "-"
                << end
                << std::endl;
        }


        // Destroy iterator
        hts_itr_destroy(itr);
    }


    // ============================================================
    // 9. Determine whether any genotype was successfully read
    // ============================================================
    bool success = !g.empty();


    if (success) {

        std::cout
            << "Genotype read in successfully. "
            << "Number of variants: "
            << g.size()
            << "; number of samples: "
            << g[0].size()
            << std::endl;

    } else {

        std::cerr
            << "Warning: no genotype records were successfully "
            << "read for the requested region(s)."
            << std::endl;
    }


    // ============================================================
    // 10. Clean up
    // ============================================================
    bcf_destroy(record);
    hts_idx_destroy(idx);
    bcf_hdr_destroy(hdr);
    bcf_close(vcf);


    // ============================================================
    // 11. Return status
    // ============================================================
    return success;
}

void read_in_genotype_trans(const std::string vcf_file, const std::vector<string> &variant_id, 
                      const std::vector<std::string> common_sample, const std::string chrom) {
    std::cout << "Start reading in Genotype!" << std::endl;
    htsFile *vcf = bcf_open(vcf_file.c_str(), "r");
    if (!vcf) {
        std::cerr << "Error opening VCF file." << std::endl;
        return;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(vcf);
    if (!hdr) {
        std::cerr << "Error reading VCF header." << std::endl;
        bcf_close(vcf);
        return;
    }

    // No index needed: we scan the entire file sequentially

    bcf1_t *record = bcf_init();
    if (!record) {
        std::cerr << "Error initializing VCF record." << std::endl;
        bcf_hdr_destroy(hdr);
        bcf_close(vcf);
        return;
    }

    std::unordered_map<std::string, int> sample_idx_map;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) {
        sample_idx_map[hdr->samples[i]] = i;
    }

    // Convert variant_id to a set for fast lookup
    std::unordered_set<std::string> variant_id_set(variant_id.begin(), variant_id.end());

    while (bcf_read(vcf, hdr, record) == 0) {
        bcf_unpack(record, BCF_UN_ALL);
        if (record->d.id && variant_id_set.find(record->d.id) != variant_id_set.end()) {

            int record_pos = record->pos + 1;
            string ref = record->d.allele[0];
            string alt = record->d.allele[1];
            string snp = chrom + ":" + to_string(record_pos) + ":" + ref + ":" + alt;
            chr_pos.push_back(snp);

            Eigen::VectorXd tmp_g = Eigen::VectorXd::Zero(common_sample.size());

            // Extract genotype array once per record, not once per sample
            int *gt_arr = nullptr, n_gts = 0;
            int gt_call_result = bcf_get_genotypes(hdr, record, &gt_arr, &n_gts);

            if (gt_call_result >= 0) {
                for (int c = 0; c < common_sample.size(); c++) {
                    string sample = common_sample[c];
                    std::string sample_substr = sample.substr(0, sample.find(":"));
                    if (sample_idx_map.find(sample_substr) != sample_idx_map.end()) {
                        int sample_index = sample_idx_map[sample_substr];
                        int dosage = 0;

                        for (int j = sample_index * 2; j < sample_index * 2 + 2 && j < n_gts; ++j) {
                            if (gt_arr[j] != bcf_gt_missing) {
                                dosage += bcf_gt_allele(gt_arr[j]);
                            }
                        }
                        tmp_g[c] = dosage;
                    }
                }
            } else {
                std::cerr << "Warning: failed to extract genotypes for record " << record->d.id
                           << " at " << chrom << ":" << record_pos << std::endl;
            }

            if (gt_arr) {
                free(gt_arr);
            }

            g.push_back(tmp_g);
        }
    }

    if (!g.empty()) {
        cout << g[0].size() << "\t" << g.size() << "\t" << g[0][1] << "\t"
             << (g.size() > 1 ? g[1][1] : 0) << "\t" << chr_pos.size() << endl;
    } else {
        std::cerr << "Warning: no matching variants found in " << vcf_file << std::endl;
    }
    cout << "Genotype read in" << endl;

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


  vector<double> dispersion_estimate(int times, Eigen::MatrixXd covariate_adjusted_geno, Eigen::VectorXd residuals, Eigen::VectorXd pi, Eigen::SparseMatrix<double> Identity){
    const double maf = 0.5;
    random_device rd;
    mt19937 gen(rd());
    binomial_distribution<> geno_dist(2, maf);
    double score_vector, info_matrix;
    Eigen::VectorXd genotype, g;
    Eigen::VectorXd orig = Eigen::VectorXd::NullaryExpr(Identity.cols(), [&]() { return geno_dist(gen); });
    double observed = 0;
    double expected = 0;

    vector<double> results;

    for(int i = 0; i < times; i++){
        for(int j = 0; j < 100; j++){
            static thread_local std::mt19937 gen(std::random_device{}());
            std::shuffle(orig.data(), orig.data() + orig.size(), gen);
            g = orig;
            genotype = Identity * g;
            // Adjust genotype
            g = genotype - covariate_adjusted_geno * genotype;

            // Compute score vector and info matrix
            score_vector = (residuals.array() * g.array()).sum();
            info_matrix = (g.array().square() * total_.array() * pi.array() * (1.0 - pi.array())).sum();

            observed = observed + score_vector * score_vector;
            expected = expected + info_matrix;
        }
    }

    double variance_ratio = std::sqrt(observed / expected);
    results.push_back(variance_ratio);

    if (convergence_status_ == false) {
        results.push_back(-1);
        return results;
    }

    return results;
}

void compute(string output_path, string site, Eigen::SparseMatrix<double> Identity, int times) {  
    // Write in middle terms for pvalue computation: residuals, pi, total, y, dispersion
    string file = output_path + "/" + site + ".middle";
    ofstream fout(file);

    Eigen::MatrixXd X_T_W_X = (X.transpose() * W_) * X;
    cout << "rows for X_T_W_X " << X_T_W_X.rows() << "\t" << X_T_W_X.cols() << endl;
    Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
    Eigen::MatrixXd X_T_W = X.transpose() * W_;
    Eigen::MatrixXd covariate_adjusted_geno = X * (X_T_W_X_inv * X_T_W);

    Eigen::VectorXd pi = sigmoid(X * beta_hat + u_hat_);
    Eigen::VectorXd mu = pi.array() * total_.array();
    Eigen::VectorXd residuals = y_.array() - mu.array();

    cout << "Start variance corrector factor estimation" << endl;

    vector<double> dispersion = dispersion_estimate(times, covariate_adjusted_geno, residuals, pi, Identity);

    string status_label = convergence_status_ ? "Converged" : "Fixed";
    fout << status_label << "\t" << site << "\t" << dispersion[0] << "\t" << tau_ << endl;

    fout << "residuals" << "\t";
    for (int i = 0; i < residuals.size(); i++) {
        fout << residuals[i] << "\t";
    }
    fout << endl;

    fout << "pi" << "\t";
    for (int i = 0; i < pi.size(); i++) {
        fout << pi[i] << "\t";
    }
    fout << endl;

    fout << "total" << "\t";
    for (int i = 0; i < total_.size(); i++) {
        fout << total_[i] << "\t";
    }
    fout << endl;

    fout << "y" << "\t";
    for (int i = 0; i < y_.size(); i++) {
        fout << y_[i] << "\t";
    }
    fout << endl;

    fout.close();
}

    void pvalue_beta_sd_compute(string output_file, string site, double dispersion, double threshold, Eigen::VectorXd residuals,Eigen::VectorXd pi,Eigen::VectorXd y,Eigen::VectorXd total){
        if(chr_pos.size()==0){
            cout<<"No genotype here"<<endl;
            return;
        }
        cout<<"Start pvalue computing!"<<endl;
        ofstream fout(output_file);
        Eigen::VectorXd tmp_t_p_1_p = total_.array() * pi.array() * (1.0 - pi.array());
        //covariate_adjusted_genotype
        Eigen::VectorXd W_vec = (total_.array()*(pi.array()*(1-pi.array()))).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
        for(int i=0;i<chr_pos.size();i++){
            const auto &g_vec = g[i];
            Eigen::VectorXd geno = g_vec - covariate_adjusted_geno * g_vec;
            // Compute score vector and info matrix
            double score_vector = (residuals.array() * geno.array()).sum();
            double info_matrix = (geno.array().square() * tmp_t_p_1_p.array() ).sum();

            // Test statistic
            double test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
            double p_value = gsl_cdf_chisq_Q(test_statistics*test_statistics,1);

    // 3. Compute effect size
            double effect = score_vector/(info_matrix*dispersion*dispersion);
    // 5. Compute standard error
            double standard_error = 1 / (std::sqrt(info_matrix) * dispersion);
            if(p_value<threshold){
                fout<<site<<"\t"<<chr_pos[i]<<"\t"<<p_value<<"\t"<<effect<<"\t"<<standard_error<<endl;}
        }
            fout.close();
    }

    void DS(string output_file, string site, double dispersion, vector<int> group, Eigen::VectorXd residuals,Eigen::VectorXd pi,Eigen::VectorXd y,Eigen::VectorXd total){
        ofstream fout(output_file);
        Eigen::VectorXd tmp_t_p_1_p = total.array() * pi.array() * (1.0 - pi.array());
        Eigen::VectorXd group_double(group.size());
        for(size_t i=0;i<group.size();i++){
            group_double[i]=static_cast<double>(group[i]);
        }
        Eigen::VectorXd W_vec = (total.array()*(pi.array()*(1-pi.array()))).matrix();
        Eigen::MatrixXd W_diag = W_vec.asDiagonal();
        W_ = W_diag.sparseView();
        Eigen::MatrixXd X_T_W_X = (X.transpose()*W_) * X;
        cout<<"rows for X_T_W_X"<<X_T_W_X.rows()<<"\t"<<X_T_W_X.rows()<<endl;
        Eigen::MatrixXd X_T_W_X_inv = X_T_W_X.inverse();
        Eigen::MatrixXd X_T_W = X.transpose() * W_;
        Eigen::MatrixXd covariate_adjusted_geno = X*(X_T_W_X_inv*X_T_W);
    // Copy values from std::vector to Eigen::VectorXd
        Eigen::VectorXd geno_double = group_double - covariate_adjusted_geno * group_double;
            // Compute score vector and info matrix
        double score_vector = (residuals.array() * geno_double.array()).sum();
        double info_matrix = (geno_double.array().square() * tmp_t_p_1_p.array() ).sum();

            // Test statistic
        double test_statistics = score_vector / (std::sqrt(info_matrix)*dispersion);
        double p_value = gsl_cdf_chisq_Q(test_statistics*test_statistics,1);
        double effect = score_vector/(info_matrix*dispersion*dispersion);
        double standard_error = 1 / (std::sqrt(info_matrix) * dispersion);

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

        // Handle different data types, using the correct pointer width for each
        if (info->type == BCF_BT_INT8) {
            int8_t* values = (int8_t*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                if (values[j] == bcf_int8_missing) {
                    infoString << ".";
                } else {
                    infoString << static_cast<int>(values[j]);
                }
            }
        } else if (info->type == BCF_BT_INT16) {
            int16_t* values = (int16_t*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                if (values[j] == bcf_int16_missing) {
                    infoString << ".";
                } else {
                    infoString << values[j];
                }
            }
        } else if (info->type == BCF_BT_INT32) {
            int32_t* values = (int32_t*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                if (values[j] == bcf_int32_missing) {
                    infoString << ".";
                } else {
                    infoString << values[j];
                }
            }
        } else if (info->type == BCF_BT_FLOAT) {
            float* values = (float*)info->vptr;
            for (int j = 0; j < info->len; ++j) {
                if (j > 0) infoString << ",";
                if (bcf_float_is_missing(values[j])) {
                    infoString << ".";
                } else {
                    infoString << values[j];
                }
            }
        } else if (info->type == BCF_BT_CHAR) {
            infoString << std::string((char*)info->vptr, info->len);
        }
    }

    return infoString.str();
}


private:
    Eigen::VectorXd& eta;
    Eigen::VectorXd& u_hat_;
    const Eigen::MatrixXd& X;
    Eigen::VectorXd& beta_hat;
    const Eigen::SparseMatrix<double> K;
    const Eigen::SparseMatrix<double> I;     
    Eigen::SparseMatrix<double> W_;
    Eigen::MatrixXd& sigma_inv_X;
    const int max_iteration_;
    Eigen::SparseMatrix<double> sigma;
    const double tol_;
    vector<Eigen::VectorXd> g;
    vector<string> chr_pos;
    vector<double> MAF;
    vector<double> distance;
    Eigen::VectorXd& y_;
    Eigen::VectorXd& total_;
    double tau_;                             
    double tau_o_;                             
    bool convergence_status_;
    Eigen::VectorXd u_hat_g_;                  
    Eigen::VectorXd u_hat_o_;                  
    vector<std::pair<double,double>> tau_trace_; 
};

int model(int argc, char *argv[]);
int QTL_parse_options(int argc, char *argv[]);
int DS_parse_options(int argc, char *argv[]);
int QTL(int argc, char *argv[]);
int trans_QTL(int argc, char *argv[]);
int DS(int argc, char *argv[]);


#endif

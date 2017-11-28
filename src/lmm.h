// linear mixed model via RcppEigen

// This code was developed following study of Nick Furlotte's pylmm code
// (https://github.com/nickFurlotte/pylmm).

#ifndef LMM_H
#define LMM_H

#include <math.h>
#include <RcppEigen.h>

struct lmm_fit {
    double hsq;
    Eigen::VectorXd beta;
    double sigmasq;
    double loglik;
    double rss;
    double logdetXSX;
};

struct calcLL_args {
    Eigen::VectorXd Kva;
    Eigen::VectorXd y;
    Eigen::MatrixXd X;
    bool reml;
    double logdetXpX;
};

struct eigenrot {
    Eigen::VectorXd Kva;
    Eigen::MatrixXd Kve;
    Eigen::MatrixXd y;
    Eigen::MatrixXd X;
};

// eigen decomposition
//    returns eigenvalues and transposed eigenvectors
std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigen_decomp(const Eigen::MatrixXd& A);

// eigen decomposition
//    returns list with eigenvalues and transposed eigenvectors
Rcpp::List Rcpp_eigen_decomp(const Rcpp::NumericMatrix &A);

// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(const Eigen::MatrixXd& K,
                               const Eigen::MatrixXd& y,
                               const Eigen::MatrixXd& X);

// eigen + rotation
Rcpp::List Rcpp_eigen_rotation(const Rcpp::NumericMatrix& K,
                               const Rcpp::NumericMatrix& y,
                               const Rcpp::NumericMatrix& X);

// calculate log det X'X
double calc_logdetXpX(const Eigen::MatrixXd& X);

// calculate log det X'X (version called from R)
double Rcpp_calc_logdetXpX(const Rcpp::NumericMatrix& X);

// getMLsoln
// for fixed value of hsq, calculate MLEs of beta and sigmasq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = whether you'll be using REML (so need to calculate log det XSX)
struct lmm_fit getMLsoln(const double hsq,
                         const Eigen::VectorXd& Kva,
                         const Eigen::VectorXd& y,
                         const Eigen::MatrixXd& X,
                         const bool reml);

// calcLL
// calculate log likelihood for fixed value of hsq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// logdetXpX = log det X'X; if NA, it's calculated
struct lmm_fit calcLL(const double hsq,
                      const Eigen::VectorXd& Kva,
                      const Eigen::VectorXd& y,
                      const Eigen::MatrixXd& X,
                      const bool reml,
                      const double logdetXpX);

// calculate log likelihood for fixed value of hsq
// This version called from R, and just returns the log likelihood
double Rcpp_calcLL(const double hsq,
                   const Rcpp::NumericVector& Kva,
                   const Rcpp::NumericVector& y,
                   const Rcpp::NumericMatrix& X,
                   const bool reml,
                   const double logdetXpX);

// just the negative log likelihood, for the optimization
double negLL(const double x, struct calcLL_args *args);

// fitLMM
// Optimize log liklihood over hsq
//
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// check_boundary = if true, explicity check 0.0 and 1.0 boundaries
// logdetXpX = log det X'X; if NA, it's calculated
// tol   = tolerance for convergence
struct lmm_fit fitLMM(const Eigen::VectorXd& Kva,
                      const Eigen::VectorXd& y,
                      const Eigen::MatrixXd& X,
                      const bool reml,
                      const bool check_boundary,
                      const double logdetXpX,
                      const double tol);

// fitLMM (version called from R)
Rcpp::List Rcpp_fitLMM(const Rcpp::NumericVector& Kva,
                       const Rcpp::NumericVector& y,
                       const Rcpp::NumericMatrix& X,
                       const bool reml,
                       const bool check_boundary,
                       const double logdetXpX,
                       const double tol);

// fitLMM with matrix of phenotypes (looping over phenotype columns)
Rcpp::List Rcpp_fitLMM_mat(const Rcpp::NumericVector& Kva, const Rcpp::NumericMatrix& Y,
                           const Rcpp::NumericMatrix& X,
                           const bool reml, const bool check_boundary,
                           const double logdetXpX, const double tol);

#endif // LMM_H

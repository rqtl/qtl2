// linear mixed model via RcppEigen

// This code was developed following study of Nick Furlotte's pylmm code
// (https://github.com/nickFurlotte/pylmm).

#ifndef LMM_H
#define LMM_H

struct lmm_fit {
    double hsq;
    VectorXd beta;
    double sigmasq;
    double loglik;
    double rss;
    double logdetXSX;
};

struct calcLL_args {
    VectorXd Kva;
    VectorXd y;
    MatrixXd X;
    bool reml;
    double logdetXpX;
};

struct eigenrot {
    VectorXd Kva;
    MatrixXd Kve;
    MatrixXd y;
    MatrixXd X;
};

// eigen decomposition
//    returns eigenvalues and transposed eigenvectors
std::pair<VectorXd, MatrixXd> eigen_decomp(const MatrixXd& A);

// eigen decomposition
//    returns list with eigenvalues and transposed eigenvectors
List Rcpp_eigen_decomp(const NumericMatrix &A);

// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(const MatrixXd& K, const MatrixXd& y,
                               const MatrixXd& X);

// eigen + rotation
List Rcpp_eigen_rotation(const NumericMatrix& K, const NumericMatrix& y,
                         const NumericMatrix& X);

// calculate log det X'X
double calc_logdetXpX(const MatrixXd& X);

// getMLsoln
// for fixed value of hsq, calculate MLEs of beta and sigmasq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = whether you'll be using REML (so need to calculate log det XSX)
struct lmm_fit getMLsoln(const double hsq, const VectorXd& Kva, const VectorXd& y,
                         const MatrixXd& X, const bool reml);

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
struct lmm_fit calcLL(const double hsq, const VectorXd& Kva, const VectorXd& y,
                const MatrixXd& X, const bool reml, const double logdetXpX);

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
struct lmm_fit fitLMM(const VectorXd& Kva, const VectorXd& y, const MatrixXd& X,
                      const bool reml, const bool check_boundary,
                      const double logdetXpX, const double tol);

// fitLMM (version called from R)
List Rcpp_fitLMM(const NumericVector& Kva, const NumericVector& y, const NumericMatrix& X,
                 const bool reml, const bool check_boundary,
                 const double logdetXpX, const double tol);

#endif // LMM_H

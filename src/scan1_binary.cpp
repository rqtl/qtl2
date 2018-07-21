// genome scan by logistic regression

#include "scan1_binary.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "binreg.h"
#include "binreg_weighted.h"
#include "matrix.h"


// Scan a single chromosome with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed, all must have values in [0,1])
// addcovar  = additive covariates (an intercept, at least)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr(const NumericVector& genoprobs,
                                 const NumericMatrix& pheno,
                                 const NumericMatrix& addcovar,
                                 const int maxit=100,
                                 const double tol=1e-6,
                                 const double qr_tol=1e-12,
                                 const double nu_max=30.0)
{
    const int n_ind = pheno.rows();
    const int n_phe = pheno.cols();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_gen = d[1];
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    const int n_add = addcovar.cols();
    const int g_size = n_ind * n_gen;

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen+n_add);

    if(n_add > 0) // paste in covariates, if present
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + g_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += g_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for this pos into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + g_size, X.begin());

        for(int phe=0; phe<n_phe; phe++) {
            // calc rss and paste into ith column of result
            result(phe,pos) = calc_ll_binreg(X, pheno(_,phe), maxit, tol, qr_tol, nu_max);
        }
    }

    return result;
}

// Scan a single chromosome with additive covariates and weights
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed, values should be in [0,1])
// addcovar  = additive covariates (an intercept, at least)
// weights   = vector of weights
//
// output    = matrix of (weighted) residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr_weighted(const NumericVector& genoprobs,
                                          const NumericMatrix& pheno,
                                          const NumericMatrix& addcovar,
                                          const NumericVector& weights,
                                          const int maxit=100,
                                          const double tol=1e-6,
                                          const double qr_tol=1e-12,
                                          const double nu_max=30.0)
{
    const int n_ind = pheno.rows();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != weights.size())
        throw std::range_error("nrow(pheno) != length(weights)");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_add = addcovar.cols();
    const int g_size = n_ind * n_gen;
    const int n_phe = pheno.cols();

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen+n_add);

    if(n_add > 0) // paste in covariates, if present
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + g_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += g_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for this pos into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + g_size, X.begin());

        for(int phe=0; phe<n_phe; phe++) {
            // calc rss and paste into ith column of result
            result(phe,pos) = calc_ll_binreg_weighted(X, pheno(_,phe), weights, maxit, tol, qr_tol, nu_max);
        }
    }

    return result;
}

// Scan a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed; values in [0,1])
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr_intcovar_highmem(const NumericVector& genoprobs,
                                                  const NumericMatrix& pheno,
                                                  const NumericMatrix& addcovar,
                                                  const NumericMatrix& intcovar,
                                                  const int maxit=100,
                                                  const double tol=1e-6,
                                                  const double qr_tol=1e-12)
{
    const int n_ind = pheno.rows();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");

    // expand genotype probabilities to include geno x interactive covariate
    NumericVector genoprobs_rev = expand_genoprobs_intcovar(genoprobs, intcovar);

    return scan_binary_onechr(genoprobs_rev, pheno, addcovar, maxit, tol, qr_tol);
}

// Scan a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr_intcovar_weighted_highmem(const NumericVector& genoprobs,
                                                           const NumericMatrix& pheno,
                                                           const NumericMatrix& addcovar,
                                                           const NumericMatrix& intcovar,
                                                           const NumericVector& weights,
                                                           const int maxit=100,
                                                           const double tol=1e-6,
                                                           const double qr_tol=1e-12)
{
    const int n_ind = pheno.rows();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_ind != weights.size())
        throw std::range_error("nrow(pheno) != length(weights)");

    // expand genotype probabilities to include geno x interactive covariate
    NumericVector genoprobs_rev = expand_genoprobs_intcovar(genoprobs, intcovar);

    // genotype can
    return scan_binary_onechr_weighted(genoprobs_rev, pheno, addcovar, weights, maxit, tol, qr_tol);
}


// Scan a single chromosome with interactive covariates
// this version uses less memory but will be slower
// (since we need to work with each position, one at a time)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr_intcovar_lowmem(const NumericVector& genoprobs,
                                                 const NumericMatrix& pheno,
                                                 const NumericMatrix& addcovar,
                                                 const NumericMatrix& intcovar,
                                                 const int maxit=100,
                                                 const double tol=1e-6,
                                                 const double qr_tol=1e-12,
                                                 const double nu_max=30.0)
{
    const int n_ind = pheno.rows();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_phe = pheno.cols();
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");

    NumericMatrix result(n_phe, n_pos);

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, true);

        for(int phe=0; phe<n_phe; phe++) {
            // do regression
            result(phe,pos) = calc_ll_binreg(X, pheno(_,phe), maxit, tol, qr_tol, nu_max);
        }
    }

    return result;
}

// Scan a single chromosome with interactive covariates
// this version uses less memory but will be slower
// (since we need to work with each position, one at a time)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_binary_onechr_intcovar_weighted_lowmem(const NumericVector& genoprobs,
                                                          const NumericMatrix& pheno,
                                                          const NumericMatrix& addcovar,
                                                          const NumericMatrix& intcovar,
                                                          const NumericVector& weights,
                                                          const int maxit=100,
                                                          const double tol=1e-6,
                                                          const double qr_tol=1e-12,
                                                          const double nu_max=30.0)
{
    const int n_ind = pheno.rows();
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs should be a 3d array but has no dim attribute");
    const Dimension d = genoprobs.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprobs should be a 3d array");
    const int n_pos = d[2];
    const int n_phe = pheno.cols();
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");

    NumericMatrix result(n_phe, n_pos);

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, true);

        for(int phe=0; phe<n_phe; phe++) {
            // do regression
            result(phe,pos) = calc_ll_binreg_weighted(X, pheno(_,phe), weights, maxit, tol, qr_tol, nu_max);
        }
    }

    return result;
}

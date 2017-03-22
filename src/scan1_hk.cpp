// genome scan by Haley-Knott regression

#include "scan1_hk.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "matrix.h"

// Scan a single chromosome with no additive covariates (not even intercept)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_nocovar(const NumericVector& genoprobs, const NumericMatrix& pheno,
                                     const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const int n_phe = pheno.cols();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int x_size = n_ind * n_gen;
    if(d[0] != n_ind)
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen);

    for(int i=0, offset=0; i<n_pos; i++, offset += x_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // calc rss and paste into ith column of result
        result(_,i) = calc_rss_linreg(X, pheno, tol);
    }

    return result;
}


// Scan a single chromosome with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr(const NumericVector& genoprobs, const NumericMatrix& pheno,
                             const NumericMatrix& addcovar, const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");

    NumericVector genoprobs_resid = calc_resid_linreg_3d(addcovar, genoprobs, tol);
    NumericMatrix pheno_resid = calc_resid_linreg(addcovar, pheno, tol);

    return scan_hk_onechr_nocovar(genoprobs_resid, pheno_resid, tol);
}

// Scan a single chromosome with additive covariates and weights
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of (weighted) residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_weighted(const NumericVector& genoprobs, const NumericMatrix& pheno,
                                      const NumericMatrix& addcovar, const NumericVector& weights,
                                      const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != weights.size())
        throw std::range_error("nrow(pheno) != length(weights)");

    // multiply everything by the (square root) of the weights
    // (weights should ALREADY be the square-root of the real weights)
    NumericMatrix addcovar_wt = weighted_matrix(addcovar, weights);
    NumericMatrix pheno_wt = weighted_matrix(pheno, weights);
    NumericVector genoprobs_wt = weighted_3darray(genoprobs, weights);

    // now regress out the additive covariates
    genoprobs_wt = calc_resid_linreg_3d(addcovar_wt, genoprobs_wt, tol);
    pheno_wt = calc_resid_linreg(addcovar_wt, pheno_wt, tol);

    // now the scan
    return scan_hk_onechr_nocovar(genoprobs_wt, pheno_wt, tol);
}

// Scan a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
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
NumericMatrix scan_hk_onechr_intcovar_highmem(const NumericVector& genoprobs,
                                              const NumericMatrix& pheno,
                                              const NumericMatrix& addcovar,
                                              const NumericMatrix& intcovar,
                                              const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");

    // expand genotype probabilities to include geno x interactive covariate
    NumericVector genoprobs_rev = expand_genoprobs_intcovar(genoprobs, intcovar);

    // regress out the additive covariates
    genoprobs_rev = calc_resid_linreg_3d(addcovar, genoprobs_rev, tol);
    NumericMatrix pheno_rev = calc_resid_linreg(addcovar, pheno, tol);

    // genotype can
    return scan_hk_onechr_nocovar(genoprobs_rev, pheno_rev, tol);
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
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_intcovar_weighted_highmem(const NumericVector& genoprobs,
                                                       const NumericMatrix& pheno,
                                                       const NumericMatrix& addcovar,
                                                       const NumericMatrix& intcovar,
                                                       const NumericVector& weights,
                                                       const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
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

    // multiply everything by the (square root) of the weights
    // (weights should ALREADY be the square-root of the real weights)
    NumericMatrix addcovar_rev = weighted_matrix(addcovar, weights);
    NumericMatrix pheno_rev = weighted_matrix(pheno, weights);
    genoprobs_rev = weighted_3darray(genoprobs_rev, weights);

    // regress out the additive covariates
    genoprobs_rev = calc_resid_linreg_3d(addcovar_rev, genoprobs_rev, tol);
    pheno_rev = calc_resid_linreg(addcovar_rev, pheno_rev, tol);

    // genotype can
    return scan_hk_onechr_nocovar(genoprobs_rev, pheno_rev, tol);
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
NumericMatrix scan_hk_onechr_intcovar_lowmem(const NumericVector& genoprobs,
                                             const NumericMatrix& pheno,
                                             const NumericMatrix& addcovar,
                                             const NumericMatrix& intcovar,
                                             const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
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

        // do regression
        result(_,pos) = calc_rss_linreg(X, pheno, tol);
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
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_intcovar_weighted_lowmem(const NumericVector& genoprobs,
                                                      const NumericMatrix& pheno,
                                                      const NumericMatrix& addcovar,
                                                      const NumericMatrix& intcovar,
                                                      const NumericVector& weights,
                                                      const double tol=1e-12)
{
    const int n_ind = pheno.rows();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_phe = pheno.cols();
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");

    NumericMatrix result(n_phe, n_pos);

    NumericMatrix pheno_rev = weighted_matrix(pheno, weights);

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, true);
        X = weighted_matrix(X, weights);

        // do regression
        result(_,pos) = calc_rss_linreg(X, pheno_rev, tol);
    }

    return result;
}

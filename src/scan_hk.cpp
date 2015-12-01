// genome scan by Haley-Knott regression

#include "scan_hk.h"
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
    const unsigned int n_ind = pheno.rows();
    const unsigned int n_phe = pheno.cols();
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_pos = d[2];
    const unsigned int n_gen = d[1];
    const unsigned int x_size = n_ind * n_gen;
    if(d[0] != n_ind)
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen);

    for(unsigned int i=0, offset=0; i<n_pos; i++, offset += x_size) {
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
    const unsigned int n_ind = pheno.rows();
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
    const unsigned int n_ind = pheno.rows();
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
    const unsigned int n_ind = pheno.rows();
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

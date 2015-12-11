// scan chromosome with linear mixed model

// [[Rcpp::depends(RcppEigen)]]

#include "scan1_lmm.h"
#include <RcppEigen.h>
#include <math.h>
#include "lmm.h"
#include "scan1_hk.h"
#include "matrix.h"
#include "linreg.h"

using namespace Rcpp;

// LMM scan of a single chromosome with additive covariates and weights
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// eigenvec  = matrix of transposed eigenvectors of variance matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = vector of log likelihood values
//
// [[Rcpp::export]]
NumericVector scan_lmm_onechr(const NumericVector& genoprobs, const NumericMatrix& pheno,
                              const NumericMatrix& addcovar, const NumericMatrix& eigenvec,
                              const NumericVector& weights, const bool reml=true,
                              const double tol=1e-12)
{
    const unsigned int n_ind = pheno.rows();
    if(pheno.cols() != 1)
        throw std::range_error("ncol(pheno) != 1");
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_gen = d[1];
    const unsigned int n_pos = d[2];
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int p = n_gen + n_addcovar;
    if(n_ind != d[0])
        throw std::range_error("ncol(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("ncol(pheno) != nrow(addcovar)");
    if(n_ind != weights.size())
        throw std::range_error("ncol(pheno) != length(weights)");
    if(n_ind != eigenvec.rows())
        throw std::range_error("ncol(pheno) != nrow(eigenvec)");
    if(n_ind != eigenvec.cols())
        throw std::range_error("ncol(pheno) != ncol(eigenvec)");

    // pre-multiply everything by the eigenvectors
    NumericVector genoprobs_copy(clone(genoprobs)); // FIX_ME: would be better not to copy
    NumericVector genoprobs_rev = matrix_x_3darray(eigenvec, genoprobs_copy);
    NumericMatrix addcovar_rev = matrix_x_matrix(eigenvec, addcovar);
    NumericMatrix pheno_rev = matrix_x_matrix(eigenvec, pheno);

    // log determinant of X'X matrices, if necessary
    NumericVector logdetXpX(n_pos);
    if(reml) logdetXpX = calc_logdetXpX_many(genoprobs_rev, addcovar_rev);

    // multiply everything by the (square root) of the weights
    // (weights should ALREADY be the square-root of the real weights)
    addcovar_rev = weighted_matrix(addcovar_rev, weights);
    pheno_rev = weighted_matrix(pheno_rev, weights);
    genoprobs_rev = weighted_3darray(genoprobs_rev, weights);

    // log determinant of X'WX matrices, if necessary
    NumericVector logdetXSX(n_pos);
    if(reml) logdetXSX = calc_logdetXpX_many(genoprobs_rev, addcovar_rev);

    // now regress out the additive covariates
    genoprobs_rev = calc_resid_linreg_3d(addcovar_rev, genoprobs_rev, tol);
    pheno_rev = calc_resid_linreg(addcovar_rev, pheno_rev, tol);

    // now the scan, return RSS
    NumericMatrix rss = scan_hk_onechr_nocovar(genoprobs_rev, pheno_rev, tol);

    // 0.5*sum(log(weights)) [since these are sqrt(weights)]
    double sum_logweights = sum(log(weights));

    NumericVector result(n_pos);
    for(unsigned int pos=0; pos<n_pos; pos++) {
        result[pos] = -(double)n_ind/2.0*log(rss[pos]) + sum_logweights;
        if(reml)
            result[pos] += 0.5*((double)p*log(2.0 * M_PI * rss[pos]/(double)(n_ind - p)) + logdetXpX[pos] - logdetXSX[pos]);
    }

    return result;
}

// REML scan of a single chromosome with interactive covariates
// this version should be fast but requires more memory
// (since we first expand the genotype probabilities to probs x intcovar)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix with one column of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = matrix of transposed eigenvectors of variance matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = vector of log likelihood values
//
// [[Rcpp::export]]
NumericVector scan_lmm_onechr_intcovar_highmem(const NumericVector& genoprobs,
                                               const NumericMatrix& pheno,
                                               const NumericMatrix& addcovar,
                                               const NumericMatrix& intcovar,
                                               const NumericMatrix& eigenvec,
                                               const NumericVector& weights,
                                               const bool reml=true,
                                               const double tol=1e-12)
{
    const unsigned int n_ind = pheno.rows();
    if(pheno.cols() != 1)
        throw std::range_error("ncol(pheno) != 1");
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_gen = d[1];
    const unsigned int n_pos = d[2];
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int n_intcovar = intcovar.cols();
    const unsigned int p = n_gen*(n_intcovar+1) + n_addcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_ind != weights.size())
        throw std::range_error("nrow(pheno) != length(weights)");
    if(n_ind != eigenvec.rows())
        throw std::range_error("ncol(pheno) != nrow(eigenvec)");
    if(n_ind != eigenvec.cols())
        throw std::range_error("ncol(pheno) != ncol(eigenvec)");

    // expand genotype probabilities to include geno x interactive covariate
    NumericVector genoprobs_rev = expand_genoprobs_intcovar(genoprobs, intcovar);

    // pre-multiply everything by eigenvectors
    genoprobs_rev = matrix_x_3darray(eigenvec, genoprobs_rev);
    NumericMatrix addcovar_rev = matrix_x_matrix(eigenvec, addcovar);
    NumericMatrix pheno_rev = matrix_x_matrix(eigenvec, pheno);

    // log determinant of X'X matrices, if necessary
    NumericVector logdetXpX(n_pos);
    if(reml) logdetXpX = calc_logdetXpX_many(genoprobs_rev, addcovar_rev);

    // multiply everything by the (square root) of the weights
    // (weights should ALREADY be the square-root of the real weights)
    addcovar_rev = weighted_matrix(addcovar_rev, weights);
    pheno_rev = weighted_matrix(pheno_rev, weights);
    genoprobs_rev = weighted_3darray(genoprobs_rev, weights);

    // log determinant of X'WX matrices, if necessary
    NumericVector logdetXSX(n_pos);
    if(reml) logdetXSX = calc_logdetXpX_many(genoprobs_rev, addcovar_rev);

    // regress out the additive covariates
    genoprobs_rev = calc_resid_linreg_3d(addcovar_rev, genoprobs_rev, tol);
    pheno_rev = calc_resid_linreg(addcovar_rev, pheno_rev, tol);

    // now the scan, return RSS
    NumericMatrix rss = scan_hk_onechr_nocovar(genoprobs_rev, pheno_rev, tol);

    // 0.5*sum(log(weights)) [since these are sqrt(weights)]
    double sum_logweights = sum(log(weights));

    // calculate log likelihood
    NumericVector result(n_pos);
    for(unsigned int pos=0; pos<n_pos; pos++) {
        result[pos] = -(double)n_ind/2.0*log(rss[pos]) + sum_logweights;
        if(reml)
            result[pos] += 0.5*((double)p*log(2.0 * M_PI * rss[pos]/(double)(n_ind - p)) + logdetXpX[pos] - logdetXSX[pos]);
    }

    return result;
}

// LMM scan of a single chromosome with interactive covariates
// this version uses less memory but will be slower
// (since we need to work with each position, one at a time)
// and this one allows weights for the individuals (the same for all phenotypes)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix with one column of numeric phenotypes
//             (no missing data allowed)
// addcovar  = additive covariates (an intercept, at least)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = matrix of transposed eigenvectors of variance matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = vector of log likelihood values
//
// [[Rcpp::export]]
NumericVector scan_lmm_onechr_intcovar_lowmem(const NumericVector& genoprobs,
                                              const NumericMatrix& pheno,
                                              const NumericMatrix& addcovar,
                                              const NumericMatrix& intcovar,
                                              const NumericMatrix& eigenvec,
                                              const NumericVector& weights,
                                              const bool reml=true,
                                              const double tol=1e-12)
{
    const unsigned int n_ind = pheno.rows();
    if(pheno.cols() != 1)
        throw std::range_error("ncol(pheno) != 1");
    const Dimension d = genoprobs.attr("dim");
    const unsigned int n_gen = d[1];
    const unsigned int n_pos = d[2];
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int n_intcovar = intcovar.cols();
    const unsigned int p = n_gen*(n_intcovar+1) + n_addcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_ind != weights.size())
        throw std::range_error("ncol(pheno) != length(weights)");
    if(n_ind != eigenvec.rows())
        throw std::range_error("ncol(pheno) != nrow(eigenvec)");
    if(n_ind != eigenvec.cols())
        throw std::range_error("ncol(pheno) != ncol(eigenvec)");

    NumericVector result(n_pos);

    // multiply pheno by eigenvectors and then by weights
    NumericMatrix pheno_rev = matrix_x_matrix(eigenvec, pheno);
    pheno_rev = weighted_matrix(pheno_rev, weights);

    // 0.5*sum(log(weights)) [since these are sqrt(weights)]
    double sum_logweights = sum(log(weights));

    for(unsigned int pos=0; pos<n_pos; pos++) {
        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos);

        // multiply by eigenvectors
        X = matrix_x_matrix(eigenvec, X);

        // log det X'X if needed
        double logdetXpX;
        if(reml) logdetXpX = Rcpp_calc_logdetXpX(X);

        // multiply by weights
        X = weighted_matrix(X, weights);

        // log det X'WX if needed
        double logdetXSX;
        if(reml) logdetXSX = Rcpp_calc_logdetXpX(X);

        // do regression
        NumericVector rss = calc_rss_linreg(X, pheno_rev, tol);
        result[pos] = -(double)n_ind/2.0*log(rss[0]) + sum_logweights;
        if(reml)
            result[pos] += 0.5*((double)p*log(2.0 * M_PI * rss[0]/(double)(n_ind - p)) + logdetXpX - logdetXSX);
    }

    return result;
}

// calculate logdetXpX many times, along positions of genotype array
NumericVector calc_logdetXpX_many(const NumericVector& genoprobs, const NumericMatrix& addcovar)
{
    if(Rf_isNull(genoprobs.attr("dim")))
        throw std::invalid_argument("genoprobs has no dimension attribute");
    const IntegerVector& dim = genoprobs.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("genoprobs should be 3-dimensional array of probabilities");
    const unsigned int n_ind = dim[0];
    const unsigned int n_gen = dim[1];
    const unsigned int n_pos = dim[2];
    const unsigned int n_addcovar = addcovar.cols();
    const unsigned int ind_by_addcovar = n_ind*n_addcovar;
    const unsigned int ind_by_gen = n_ind*n_gen;
    if(addcovar.rows() != n_ind)
        throw std::range_error("nrow(genoprobs) != nrow(addcovar)");

    NumericVector result(n_pos);
    NumericMatrix X(n_ind, n_addcovar+n_gen);
    std::copy(addcovar.begin(), addcovar.end(), X.begin());

    for(unsigned int pos=0; pos<n_pos; pos++) {
        std::copy(genoprobs.begin()+ind_by_gen*pos, genoprobs.begin()+ind_by_gen*(pos+1),
                  X.begin()+ind_by_addcovar);
        result[pos] = Rcpp_calc_logdetXpX(X);
    }

    return result;
}

// fit a single-QTL model at a single position by LMM

#include "fit1_pg.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "matrix.h"

// fit single-QTL model at a single position
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list with a bunch of stuff
//
// [[Rcpp::export]]
List fit1_pg_addcovar(const NumericMatrix& genoprobs,
                      const NumericVector& pheno,
                      const NumericMatrix& addcovar,
                      const NumericMatrix& eigenvec,
                      const NumericVector& weights,
                      const bool se=false,
                      const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const int n_gen = genoprobs.cols();
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != genoprobs.rows())
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    NumericMatrix X(n_ind, n_coef);

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;

    // copy genoprobs into matrix
    std::copy(genoprobs.begin(), genoprobs.end(), X.begin());

    // copy addcovar into matrix
    if(n_addcovar > 0)
        std::copy(addcovar.begin(), addcovar.end(), X.begin() + x_size);

    // pre-multiply by eigenvec then multiply by weights
    X = matrix_x_matrix(eigenvec, X);
    X = weighted_matrix(X, weights);

    // do regression
    List result = fit_linreg(X, pheno_rev, se, tol);

    // fix the fitted values (leave the residuals as they are)
    NumericVector fitted = result["fitted"];
    NumericVector fitted_rev = matrix_x_vector(transpose(eigenvec), fitted/weights);

    result["fitted"] = fitted_rev;

    return result;
}


// fit single-QTL model at a single position
//
// genoprobs = matrix of genotype probabilities (individuals x genotypes)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of a bunch of stuff
//
// [[Rcpp::export]]
List fit1_pg_intcovar(const NumericMatrix& genoprobs,
                      const NumericVector& pheno,
                      const NumericMatrix& addcovar,
                      const NumericMatrix& intcovar,
                      const NumericMatrix& eigenvec,
                      const NumericVector& weights,
                      const bool se=true,
                      const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const int n_weights = weights.size();
    if(n_ind != genoprobs.rows())
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;

    // form X matrix
    NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, 0, false);

    // pre-multiply by eigenvec then multiply by weights
    X = matrix_x_matrix(eigenvec, X);
    X = weighted_matrix(X, weights);

    // do regression
    List result = fit_linreg(X, pheno_rev, se, tol);

    // fix the fitted values (leave the residuals as they are)
    NumericVector fitted = result["fitted"];
    NumericVector fitted_rev = matrix_x_vector(transpose(eigenvec), fitted/weights);

    result["fitted"] = fitted_rev;

    return result;
}

// scan chromosome by LMM (to adjust for polygenic effect) just to get coefficients

#include "scan1coef_pg.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "matrix.h"

// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scancoef_pg_addcovar(const NumericVector& genoprobs,
                                   const NumericVector& pheno,
                                   const NumericMatrix& addcovar,
                                   const NumericMatrix& eigenvec,
                                   const NumericVector& weights,
                                   const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    NumericMatrix result(n_coef, n_pos);
    NumericMatrix X(n_ind, n_coef);

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;
    NumericVector genoprobs_copy(clone(genoprobs)); // FIX_ME: would be better not to copy
    NumericVector genoprobs_rev = matrix_x_3darray(eigenvec, genoprobs_copy);
    genoprobs_rev = weighted_3darray(genoprobs_rev, weights);
    NumericMatrix addcovar_rev;
    if(n_addcovar > 0)  {
        addcovar_rev = matrix_x_matrix(eigenvec, addcovar);
        addcovar_rev = weighted_matrix(addcovar_rev, weights);
    }

    // copy addcovar into matrix
    if(n_addcovar > 0)
        std::copy(addcovar_rev.begin(), addcovar_rev.end(), X.begin() + x_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs_rev.begin() + offset, genoprobs_rev.begin() + offset + x_size, X.begin());

        // do regression
        result(_,pos) = calc_coef_linreg(X, pheno_rev, tol);
    }

    return result;
}


// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// intcovar  = interactive covariates (should also be included in addcovar)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = matrix of coefficients (genotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scancoef_pg_intcovar(const NumericVector& genoprobs,
                                   const NumericVector& pheno,
                                   const NumericMatrix& addcovar,
                                   const NumericMatrix& intcovar,
                                   const NumericMatrix& eigenvec,
                                   const NumericVector& weights,
                                   const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int n_intcovar = intcovar.cols();
    const int n_coef = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    NumericMatrix result(n_coef, n_pos);

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, false);

        // pre-multiply by eigenvec then multiply by weights
        X = matrix_x_matrix(eigenvec, X);
        X = weighted_matrix(X, weights);

        // do regression
        result(_,pos) = calc_coef_linreg(X, pheno_rev, tol);
    }

    return result;
}


// Scan a single chromosome to calculate coefficients, with additive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of two matrices, of coefficients and SEs (each genotypes x positions)
//
// [[Rcpp::export]]
List scancoefSE_pg_addcovar(const NumericVector& genoprobs,
                            const NumericVector& pheno,
                            const NumericMatrix& addcovar,
                            const NumericMatrix& eigenvec,
                            const NumericVector& weights,
                            const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int x_size = n_ind * n_gen;
    const int n_coef = n_gen + n_addcovar;

    if(n_ind != d[0])
        throw std::range_error("length(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("length(pheno) != nrow(addcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    NumericMatrix coef(n_coef, n_pos);
    NumericMatrix se(n_coef, n_pos);
    NumericMatrix X(n_ind, n_coef);

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;
    NumericVector genoprobs_copy(clone(genoprobs)); // FIX_ME: would be better not to copy
    NumericVector genoprobs_rev = matrix_x_3darray(eigenvec, genoprobs_copy);
    genoprobs_rev = weighted_3darray(genoprobs_rev, weights);
    NumericMatrix addcovar_rev;
    if(n_addcovar > 0) {
        addcovar_rev = matrix_x_matrix(eigenvec, addcovar);
        addcovar_rev = weighted_matrix(addcovar_rev, weights);
    }

    // copy addcovar into matrix
    if(n_addcovar > 0)
        std::copy(addcovar_rev.begin(), addcovar_rev.end(), X.begin() + x_size);

    for(int pos=0, offset=0; pos<n_pos; pos++, offset += x_size) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs_rev.begin() + offset, genoprobs_rev.begin() + offset + x_size, X.begin());

        // do regression
        List tmp = calc_coefSE_linreg(X, pheno_rev, tol);
        NumericVector tmpcoef = tmp[0];
        NumericVector tmpse = tmp[1];
        coef(_,pos) = tmpcoef;
        se(_,pos) = tmpse;
    }

    return List::create(Named("coef") = coef,
                        Named("SE") = se);
}


// Scan a single chromosome to calculate coefficients, with interactive covariates
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = vector of numeric phenotypes (individuals x 1)
//             (no missing data allowed)
// addcovar  = additive covariates (can be null)
// eigenvec  = eigenvectors from eigen decomposition of kinship matrix
// intcovar  = interactive covariates (should also be included in addcovar)
// weights   = vector of weights (really the SQUARE ROOT of the weights)
//
// output    = list of two matrices, of coefficients and SEs (each genotypes x positions)
//
// [[Rcpp::export]]
List scancoefSE_pg_intcovar(const NumericVector& genoprobs,
                            const NumericVector& pheno,
                            const NumericMatrix& addcovar,
                            const NumericMatrix& intcovar,
                            const NumericMatrix& eigenvec,
                            const NumericVector& weights,
                            const double tol=1e-12)
{
    const int n_ind = pheno.size();
    const Dimension d = genoprobs.attr("dim");
    const int n_pos = d[2];
    const int n_gen = d[1];
    const int n_weights = weights.size();
    const int n_addcovar = addcovar.cols();
    const int n_intcovar = intcovar.cols();
    const int n_coef = n_gen + n_addcovar + (n_gen-1)*n_intcovar;
    if(n_ind != d[0])
        throw std::range_error("nrow(pheno) != nrow(genoprobs)");
    if(n_ind != addcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(addcovar)");
    if(n_ind != intcovar.rows())
        throw std::range_error("nrow(pheno) != nrow(intcovar)");
    if(n_weights != n_ind)
        throw std::range_error("length(pheno) != length(weights)");
    if(eigenvec.rows() != n_ind || eigenvec.cols() != n_ind)
        throw std::range_error("eigenvec should be square matrix with dimension length(pheno)");

    NumericMatrix coef(n_coef, n_pos);
    NumericMatrix se(n_coef, n_pos);

    // pre-multiply by eigenvectors then multiply by weights
    NumericVector pheno_rev = matrix_x_vector(eigenvec, pheno);
    pheno_rev = pheno_rev * weights;

    for(int pos=0; pos<n_pos; pos++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // form X matrix
        NumericMatrix X = formX_intcovar(genoprobs, addcovar, intcovar, pos, false);

        // pre-multiply by eigenvec then multiply by weights
        X = matrix_x_matrix(eigenvec, X);
        X = weighted_matrix(X, weights);

        // do regression
        List tmp = calc_coefSE_linreg(X, pheno_rev, tol);
        NumericVector tmpcoef = tmp[0];
        NumericVector tmpse = tmp[1];
        coef(_,pos) = tmpcoef;
        se(_,pos) = tmpse;
    }

    return List::create(Named("coef") = coef,
                        Named("SE") = se);
}

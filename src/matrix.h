// Matrix utilities
#ifndef MATRIX_H
#define MATRIX_H

#include <RcppEigen.h>

// cbind two or three matrices
Rcpp::IntegerMatrix cbind_imatrix(const Rcpp::IntegerMatrix& mat1,
                                  const Rcpp::IntegerMatrix& mat2);
Rcpp::IntegerMatrix cbind_3imatrix(const Rcpp::IntegerMatrix& mat1,
                                   const Rcpp::IntegerMatrix& mat2,
                                   const Rcpp::IntegerMatrix& mat3);
Rcpp::NumericMatrix cbind_nmatrix(const Rcpp::NumericMatrix& mat1,
                                  const Rcpp::NumericMatrix& mat2);
Rcpp::NumericMatrix cbind_3nmatrix(const Rcpp::NumericMatrix& mat1,
                                   const Rcpp::NumericMatrix& mat2,
                                   const Rcpp::NumericMatrix& mat3);

// rbind two or three matrices
Rcpp::IntegerMatrix rbind_imatrix(const Rcpp::IntegerMatrix& mat1,
                                  const Rcpp::IntegerMatrix& mat2);
Rcpp::IntegerMatrix rbind_3imatrix(const Rcpp::IntegerMatrix& mat1,
                                   const Rcpp::IntegerMatrix& mat2,
                                   const Rcpp::IntegerMatrix& mat3);
Rcpp::NumericMatrix rbind_nmatrix(const Rcpp::NumericMatrix& mat1,
                                  const Rcpp::NumericMatrix& mat2);
Rcpp::NumericMatrix rbind_3nmatrix(const Rcpp::NumericMatrix& mat1,
                                   const Rcpp::NumericMatrix& mat2,
                                   const Rcpp::NumericMatrix& mat3);

// find columns that exactly match previous columns
// returns numeric vector with -1 indicating no match to an earlier column and
//                             >0 indicating matches that earlier column
//                                (indexes starting at 1)
Rcpp::NumericVector find_matching_cols(const Rcpp::NumericMatrix& mat,
                                       const double tol);

// find set of linearly independent columns in a matrix
// returns a vector of column indices (starting at 1)
Rcpp::IntegerVector find_lin_indep_cols(const Rcpp::NumericMatrix& mat,
                                        const double tol);

// form X matrix with intcovar
Rcpp::NumericMatrix formX_intcovar(const Rcpp::NumericMatrix& probs,
                                   const Rcpp::NumericMatrix& addcovar,
                                   const Rcpp::NumericMatrix& intcovar);

// multiply each column of a matrix by a set of weights
Rcpp::NumericMatrix weighted_matrix(const Rcpp::NumericMatrix& mat,
                              const Rcpp::NumericVector& weights);

// multiply each element of a vector by the corresponding weight
Rcpp::NumericVector weighted_3darray(const Rcpp::NumericVector& array,
                                     const Rcpp::NumericVector& weights);

#endif // MATRIX_H

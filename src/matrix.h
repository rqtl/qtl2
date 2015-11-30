// Matrix utilities
#ifndef MATRIX_H
#define MATRIX_H

// cbind two or three matrices
IntegerMatrix cbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2);
IntegerMatrix cbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3);
NumericMatrix cbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2);
NumericMatrix cbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3);

// rbind two or three matrices
IntegerMatrix rbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2);
IntegerMatrix rbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3);
NumericMatrix rbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2);
NumericMatrix rbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3);

// find columns that exactly match previous columns
// returns numeric vector with -1 indicating no match to an earlier column and
//                             >0 indicating matches that earlier column
//                                (indexes starting at 1)
NumericVector find_matching_cols(const NumericMatrix& mat, const double tol);

// find set of linearly independent columns in a matrix
// returns a vector of column indices (starting at 1)
IntegerVector find_lin_indep_cols(const NumericMatrix& mat, const double tol);

#endif // MATRIX_H

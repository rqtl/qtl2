// Matrix utilities

#include <Rcpp.h>
using namespace Rcpp;

#include "matrix.h"

// cbind two matrices
// [[Rcpp::export]]
IntegerMatrix cbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2)
{
    unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    unsigned int ncol1 = mat1.cols();
    unsigned int ncol2 = mat2.cols();

    unsigned int ncol = ncol1 + ncol2;
    IntegerMatrix result(nrow, ncol);

    std::copy(mat1.begin(), mat1.end(), result.begin());
    unsigned int curindex = nrow*ncol1;
    std::copy(mat2.begin(), mat2.end(), result.begin() + curindex);

    return result;
}

// cbind three matrices
// [[Rcpp::export]]
IntegerMatrix cbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3)
{
    unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    if(mat3.rows() != nrow)
        throw std::length_error("nrow(mat3) != nrow(mat1)");
    unsigned int ncol1 = mat1.cols();
    unsigned int ncol2 = mat2.cols();
    unsigned int ncol3 = mat3.cols();

    unsigned int ncol = ncol1 + ncol2 + ncol3;
    IntegerMatrix result(nrow, ncol);

    std::copy(mat1.begin(), mat1.end(), result.begin());
    unsigned int curindex = nrow*ncol1;
    std::copy(mat2.begin(), mat2.end(), result.begin() + curindex);
    curindex += nrow*ncol2;
    std::copy(mat3.begin(), mat3.end(), result.begin() + curindex);

    return result;
}

// cbind two matrices
// [[Rcpp::export]]
NumericMatrix cbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2)
{
    unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    unsigned int ncol1 = mat1.cols();
    unsigned int ncol2 = mat2.cols();

    unsigned int ncol = ncol1 + ncol2;
    NumericMatrix result(nrow, ncol);

    std::copy(mat1.begin(), mat1.end(), result.begin());
    unsigned int curindex = nrow*ncol1;
    std::copy(mat2.begin(), mat2.end(), result.begin() + curindex);

    return result;
}

// cbind three matrices
// [[Rcpp::export]]
NumericMatrix cbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3)
{
    unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    if(mat3.rows() != nrow)
        throw std::length_error("nrow(mat3) != nrow(mat1)");
    unsigned int ncol1 = mat1.cols();
    unsigned int ncol2 = mat2.cols();
    unsigned int ncol3 = mat3.cols();

    unsigned int ncol = ncol1 + ncol2 + ncol3;
    NumericMatrix result(nrow, ncol);

    std::copy(mat1.begin(), mat1.end(), result.begin());
    unsigned int curindex = nrow*ncol1;
    std::copy(mat2.begin(), mat2.end(), result.begin() + curindex);
    curindex += nrow*ncol2;
    std::copy(mat3.begin(), mat3.end(), result.begin() + curindex);

    return result;
}

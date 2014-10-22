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

// rbind two matrices
// [[Rcpp::export]]
IntegerMatrix rbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2)
{
    unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    unsigned int nrow1 = mat1.rows();
    unsigned int nrow2 = mat2.rows();

    unsigned int nrow = nrow1 + nrow2;
    IntegerMatrix result(nrow, ncol);

    for(unsigned int col=0, offset1=0, offset2=0, offset_result=0;
        col < ncol;
        col++, offset1 += nrow1, offset2 += nrow2, offset_result += nrow) {

        std::copy(mat1.begin()+offset1, mat1.begin()+offset1+nrow1, result.begin()+offset_result);
        std::copy(mat2.begin()+offset2, mat2.begin()+offset2+nrow2, result.begin()+offset_result+nrow1);
    }

    return result;
}

// rbind three matrices
// [[Rcpp::export]]
IntegerMatrix rbind_3imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2,
                             const IntegerMatrix& mat3)
{
    unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    if(mat3.cols() != ncol)
        throw std::length_error("ncol(mat3) != ncol(mat1)");
    unsigned int nrow1 = mat1.rows();
    unsigned int nrow2 = mat2.rows();
    unsigned int nrow3 = mat3.rows();

    unsigned int nrow = nrow1 + nrow2 + nrow3;
    IntegerMatrix result(nrow, ncol);

    for(unsigned int col=0, offset1=0, offset2=0, offset3=0, offset_result=0;
        col < ncol;
        col++, offset1 += nrow1, offset2 += nrow2, offset3 += nrow3, offset_result += nrow) {

        std::copy(mat1.begin()+offset1, mat1.begin()+offset1+nrow1, result.begin()+offset_result);
        std::copy(mat2.begin()+offset2, mat2.begin()+offset2+nrow2, result.begin()+offset_result+nrow1);
        std::copy(mat3.begin()+offset3, mat3.begin()+offset3+nrow3, result.begin()+offset_result+nrow1+nrow2);
    }

    return result;
}

// rbind two matrices
// [[Rcpp::export]]
NumericMatrix rbind_nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2)
{
    unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    unsigned int nrow1 = mat1.rows();
    unsigned int nrow2 = mat2.rows();

    unsigned int nrow = nrow1 + nrow2;
    NumericMatrix result(nrow, ncol);

    for(unsigned int col=0, offset1=0, offset2=0, offset_result=0;
        col < ncol;
        col++, offset1 += nrow1, offset2 += nrow2, offset_result += nrow) {

        std::copy(mat1.begin()+offset1, mat1.begin()+offset1+nrow1, result.begin()+offset_result);
        std::copy(mat2.begin()+offset2, mat2.begin()+offset2+nrow2, result.begin()+offset_result+nrow1);
    }

    return result;
}

// rbind three matrices
// [[Rcpp::export]]
NumericMatrix rbind_3nmatrix(const NumericMatrix& mat1, const NumericMatrix& mat2,
                             const NumericMatrix& mat3)
{
    unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    if(mat3.cols() != ncol)
        throw std::length_error("ncol(mat3) != ncol(mat1)");
    unsigned int nrow1 = mat1.rows();
    unsigned int nrow2 = mat2.rows();
    unsigned int nrow3 = mat3.rows();

    unsigned int nrow = nrow1 + nrow2 + nrow3;
    NumericMatrix result(nrow, ncol);

    for(unsigned int col=0, offset1=0, offset2=0, offset3=0, offset_result=0;
        col < ncol;
        col++, offset1 += nrow1, offset2 += nrow2, offset3 += nrow3, offset_result += nrow) {

        std::copy(mat1.begin()+offset1, mat1.begin()+offset1+nrow1, result.begin()+offset_result);
        std::copy(mat2.begin()+offset2, mat2.begin()+offset2+nrow2, result.begin()+offset_result+nrow1);
        std::copy(mat3.begin()+offset3, mat3.begin()+offset3+nrow3, result.begin()+offset_result+nrow1+nrow2);
    }

    return result;
}

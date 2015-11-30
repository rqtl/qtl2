// Matrix utilities

#include <Rcpp.h>
using namespace Rcpp;

#include "matrix.h"

// cbind two matrices
// [[Rcpp::export]]
IntegerMatrix cbind_imatrix(const IntegerMatrix& mat1, const IntegerMatrix& mat2)
{
    const unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    const unsigned int ncol1 = mat1.cols();
    const unsigned int ncol2 = mat2.cols();

    const unsigned int ncol = ncol1 + ncol2;
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
    const unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    if(mat3.rows() != nrow)
        throw std::length_error("nrow(mat3) != nrow(mat1)");
    const unsigned int ncol1 = mat1.cols();
    const unsigned int ncol2 = mat2.cols();
    const unsigned int ncol3 = mat3.cols();

    const unsigned int ncol = ncol1 + ncol2 + ncol3;
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
    const unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    const unsigned int ncol1 = mat1.cols();
    const unsigned int ncol2 = mat2.cols();

    const unsigned int ncol = ncol1 + ncol2;
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
    const unsigned int nrow = mat1.rows();
    if(mat2.rows() != nrow)
        throw std::length_error("nrow(mat2) != nrow(mat1)");
    if(mat3.rows() != nrow)
        throw std::length_error("nrow(mat3) != nrow(mat1)");
    const unsigned int ncol1 = mat1.cols();
    const unsigned int ncol2 = mat2.cols();
    const unsigned int ncol3 = mat3.cols();

    const unsigned int ncol = ncol1 + ncol2 + ncol3;
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
    const unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    const unsigned int nrow1 = mat1.rows();
    const unsigned int nrow2 = mat2.rows();

    const unsigned int nrow = nrow1 + nrow2;
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
    const unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    if(mat3.cols() != ncol)
        throw std::length_error("ncol(mat3) != ncol(mat1)");
    const unsigned int nrow1 = mat1.rows();
    const unsigned int nrow2 = mat2.rows();
    const unsigned int nrow3 = mat3.rows();

    const unsigned int nrow = nrow1 + nrow2 + nrow3;
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
    const unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    const unsigned int nrow1 = mat1.rows();
    const unsigned int nrow2 = mat2.rows();

    const unsigned int nrow = nrow1 + nrow2;
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
    const unsigned int ncol = mat1.cols();
    if(mat2.cols() != ncol)
        throw std::length_error("ncol(mat2) != ncol(mat1)");
    if(mat3.cols() != ncol)
        throw std::length_error("ncol(mat3) != ncol(mat1)");
    const unsigned int nrow1 = mat1.rows();
    const unsigned int nrow2 = mat2.rows();
    const unsigned int nrow3 = mat3.rows();

    const unsigned int nrow = nrow1 + nrow2 + nrow3;
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

// find columns that exactly match previous columns
// returns numeric vector with -1 indicating no match to an earlier column and
//                             >0 indicating matches that earlier column
//                                (indexes starting at 1)
// [[Rcpp::export]]
NumericVector find_matching_cols(const NumericMatrix& mat, const double tol=1e-8)
{
    const unsigned int ncol = mat.cols();
    const unsigned int nrow = mat.rows();
    NumericVector result(ncol);

    if(ncol < 1) Rf_error("Matrix has 0 columns");

    result[0] = -1;
    if(ncol==1) return(result);

    for(unsigned int i=1; i<ncol; i++) {
        result[i] = -1;
        for(unsigned int j=0; j<i; j++) {
            double max_diff=0.0;
            for(unsigned int k=0; k<nrow; k++) {
                const bool na_i = NumericVector::is_na(mat(k,i));
                const bool na_j = NumericVector::is_na(mat(k,j));
                // if both missing, return 0.0
                // if one missing but not other, return 1.0
                // otherwise, return difference
                double d = na_i != na_j ? 1.0 : ((na_i && na_j) ? 0.0 : fabs(mat(k,i) - mat(k,j)));
                if(d > max_diff) max_diff = d;
            }
            if(max_diff < tol) {
                result[i] = j+1;
                break;
            }
        }
    }

    return(result);
}

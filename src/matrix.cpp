// Matrix utilities

// [[Rcpp::depends(RcppEigen)]]

#include "matrix.h"
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;


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
NumericVector find_matching_cols(const NumericMatrix& mat, const double tol=1e-12)
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

    return result;
}

// find set of linearly independent columns in a matrix
// returns a vector of column indices (starting at 1)
// [[Rcpp::export]]
IntegerVector find_lin_indep_cols(const NumericMatrix& mat, const double tol=1e-12)
{
    const unsigned int ncol=mat.cols();

    // QR decomp with column pivoting
    MatrixXd XX(as<Map<MatrixXd> >(mat));
    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;
    CPivQR PQR = XX;
    PQR.setThreshold(tol);

    // pivot matrix, treated as regular matrix
    Permutation Pmat = PQR.colsPermutation();
    MatrixXd PPmat(Pmat);

    // rank of input matrix
    const unsigned int rank=PQR.rank();
    IntegerVector result(rank);

    // for each column, find the row with a 1
    for(unsigned int j=0; j<rank; j++) {
        for(unsigned int i=0; i<ncol; i++) {
            if(fabs(PPmat(i,j) - 1.0) < tol) {
                result[j] = i+1;
                break;
            }
        }
    }

    return result;
}

// form X matrix with intcovar
// [[Rcpp::export]]
NumericMatrix formX_intcovar(const NumericMatrix& probs,
                             const NumericMatrix& addcovar,
                             const NumericMatrix& intcovar)
{
    const unsigned int nrow  = probs.rows();
    const unsigned int nprob = probs.cols();
    const unsigned int nadd  = addcovar.cols();
    const unsigned int nint  = intcovar.cols();

    if(addcovar.rows() != nrow)
        throw std::range_error("nrow(addcovar) != nrow(probs)");
    if(intcovar.rows() != nrow)
        throw std::range_error("nrow(intcovar) != nrow(probs)");

    NumericMatrix result(nrow,nadd+nprob+nprob*nint);
    std::copy(addcovar.begin(), addcovar.end(), result.begin());
    std::copy(probs.begin(), probs.end(), result.begin() + nrow*nadd);

    for(unsigned int i=0, rescol=nprob+nadd; i<nint; i++) {
        for(unsigned int j=0; j<nprob; j++, rescol++) {
            for(unsigned int k=0; k<nrow; k++)
                result(k,rescol) = probs(k,j)*intcovar(k,i);
        }
    }

    return result;
}

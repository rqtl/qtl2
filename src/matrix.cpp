// Matrix utilities

// [[Rcpp::depends(RcppEigen)]]

#include "matrix.h"
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

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
    const MatrixXd XX(as<Map<MatrixXd> >(mat));
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
NumericMatrix formX_intcovar(const NumericVector& probs,
                             const NumericMatrix& addcovar,
                             const NumericMatrix& intcovar,
                             const int position) // with indexes starting at 0
{

    const Dimension d = probs.attr("dim");
    const unsigned int nrow  = d[0];
    const unsigned int ngen = d[1];
    const unsigned int recsize = nrow*ngen;
    const unsigned int offset = recsize*position;
    const unsigned int nadd  = addcovar.cols();
    const unsigned int nint  = intcovar.cols();

    NumericMatrix result(nrow, nadd + ngen*(nint+1));

    if(position < 0 || position >= d[2])
        throw std::range_error("position out of range of 0 .. (n_pos-1)");
    if(addcovar.rows() != nrow)
        throw std::range_error("nrow(addcovar) != nrow(probs)");
    if(intcovar.rows() != nrow)
        throw std::range_error("nrow(intcovar) != nrow(probs)");

    std::copy(addcovar.begin(), addcovar.end(), result.begin());
    std::copy(probs.begin()+offset,
              probs.begin()+offset+recsize,
              result.begin() + nrow*nadd);

    for(unsigned int i=0, rescol=ngen+nadd; i<nint; i++) {
        for(unsigned int j=0; j<ngen; j++, rescol++) {
            for(unsigned int k=0; k<nrow; k++) {
                result(k,rescol) = probs[offset + k + j*nrow] * intcovar(k,i);
            }
        }
    }

    return result;
}


// expand genotype probabilities with intcovar
// [[Rcpp::export]]
NumericVector expand_genoprobs_intcovar(const NumericVector& probs, // 3d array ind x prob x pos
                                        const NumericMatrix& intcovar)
{
    Dimension d = probs.attr("dim");
    const unsigned int nrow  = d[0];
    const unsigned int ngen = d[1];
    const unsigned int npos = d[2];
    const unsigned int nint  = intcovar.cols();

    if(intcovar.rows() != nrow)
        throw std::range_error("nrow(intcovar) != nrow(probs)");

    const unsigned int ngen_result = d[1]*(nint+1); // no. cols in result
    const unsigned int recsize = nrow*ngen; // ind x geno rectangle
    const unsigned int recsize_result = nrow*ngen_result; // ind x geno rectangle in result

    NumericVector result(recsize_result*npos);

    for(unsigned int i=0; i<npos; i++) {
        // paste probs into first batch
        std::copy(probs.begin()+i*recsize,
                  probs.begin()+(i+1)*recsize,
                  result.begin()+i*recsize_result);
        for(unsigned int j=0; j<nint; j++) {
            for(unsigned int k=0; k<ngen; k++) {
                for(unsigned int s=0; s<nrow; s++)
                    result[i*recsize_result + (j+1)*recsize + k*nrow + s] =
                        probs[i*recsize + k*nrow + s] * intcovar(s,j);
            }
        }
    }

    // add dimension attribute
    d[1] = ngen_result;
    result.attr("dim") = d;
    rownames(result) = rownames(probs);
    return result;
}

// multiply each column of a matrix by a set of weights
// [[Rcpp::export]]
NumericMatrix weighted_matrix(const NumericMatrix& mat,
                              const NumericVector& weights)
{
    const unsigned int nrow = mat.rows();
    const unsigned int ncol = mat.cols();
    if(nrow != weights.size())
        throw std::range_error("nrow(mat) != length(weights)");

    NumericMatrix result(nrow,ncol);

    for(unsigned int j=0; j<ncol; j++)
        for(unsigned int i=0; i<nrow; i++)
            result(i,j) = mat(i,j)*weights[i];

    return result;
}

// multiply each column of a 3-dimensional array by a set of weights
// [[Rcpp::export]]
NumericVector weighted_3darray(const NumericVector& array,
                               const NumericVector& weights)
{
    const Dimension d = array.attr("dim");
    const unsigned int n = d[0];
    const unsigned int ncol = d[1]*d[2];
    if(n != weights.size())
        throw std::range_error("nrow(array) != length(weights)");

    NumericVector result(n*ncol);
    result.attr("dim") = d;

    for(unsigned int j=0, k=0; j<ncol; j++)
        for(unsigned int i=0; i<n; i++, k++)
            result[k] = array[k]*weights[i];

    return result;
}

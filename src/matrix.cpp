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
    const int ncol = mat.cols();
    const int nrow = mat.rows();
    NumericVector result(ncol);

    if(ncol < 1) Rf_error("Matrix has 0 columns");

    result[0] = -1;
    if(ncol==1) return(result);

    for(int i=1; i<ncol; i++) {
        result[i] = -1;
        for(int j=0; j<i; j++) {
            double max_diff=0.0;
            for(int k=0; k<nrow; k++) {
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
    const int ncol=mat.cols();

    // QR decomp with column pivoting
    const MatrixXd XX(as<Map<MatrixXd> >(mat));
    typedef Eigen::ColPivHouseholderQR<MatrixXd> CPivQR;
    typedef CPivQR::PermutationType Permutation;
    CPivQR PQR ( XX );
    PQR.setThreshold(tol);

    // pivot matrix, treated as regular matrix
    Permutation Pmat = PQR.colsPermutation();
    MatrixXd PPmat(Pmat);

    // rank of input matrix
    const int rank=PQR.rank();
    IntegerVector result(rank);

    // for each column, find the row with a 1
    for(int j=0; j<rank; j++) {
        for(int i=0; i<ncol; i++) {
            if(fabs(PPmat(i,j) - 1.0) < tol) {
                result[j] = i+1;
                break;
            }
        }
    }

    return result;
}

// form X matrix with intcovar
// has_intercept = true indicates that addcovar has an intercept
//                 and so probs reduced by one column
//               = false means probs has full set of columns
//
// This is maybe a bit confusing.
//
// In the has_intercept=true case, an intercept is included in the
// addcovar matrix, and probs is missing the first column.
// The matrix formed is [A P (P.I)] where A=addcovar, P=probs, I=intcovar
//
// In the has_intercept=false case, no intercept is included in the
// addcovar matrix, and the probs have all columns (so each row sums
// to 1). The matrix formed is [P A (P*.I)] where in the P*.I bit we
// drop the first column of the probs when getting interactions with
// intercovar.
//
// [[Rcpp::export]]
NumericMatrix formX_intcovar(const NumericVector& probs,
                             const NumericMatrix& addcovar,
                             const NumericMatrix& intcovar,
                             const int position, // with indexes starting at 0
                             const bool has_intercept=true)
{

    const Dimension d = probs.attr("dim");
    const int nrow  = d[0];
    const int ngen = d[1];
    const int recsize = nrow*ngen;
    const int offset = recsize*position;
    const int nadd  = addcovar.cols();
    const int nint  = intcovar.cols();

    int totcol;
    if(has_intercept) totcol = nadd + ngen*(nint+1);
    else totcol = ngen + nadd + (ngen-1)*nint;

    NumericMatrix result(nrow, totcol);

    if(position < 0 || (position != 0 && position >= d[2])) // if position == 0, don't look at d[2]
        throw std::range_error("position out of range of 0 .. (n_pos-1)");
    if(addcovar.rows() != nrow)
        throw std::range_error("nrow(addcovar) != nrow(probs)");
    if(intcovar.rows() != nrow)
        throw std::range_error("nrow(intcovar) != nrow(probs)");

    if(has_intercept) {
        std::copy(addcovar.begin(), addcovar.end(), result.begin());
        std::copy(probs.begin()+offset,
                  probs.begin()+offset+recsize,
                  result.begin() + nrow*nadd);

        for(int i=0, rescol=ngen+nadd; i<nint; i++) {
            for(int j=0; j<ngen; j++, rescol++) {
                for(int k=0; k<nrow; k++) {
                    result(k,rescol) = probs[offset + k + j*nrow] * intcovar(k,i);
                }
            }
        }
    }
    else {
        std::copy(probs.begin()+offset,
                  probs.begin()+offset+recsize,
                  result.begin());
        std::copy(addcovar.begin(), addcovar.end(), result.begin() + nrow*ngen);

        for(int i=0, rescol=ngen+nadd; i<nint; i++) {
            for(int j=1; j<ngen; j++, rescol++) {
                for(int k=0; k<nrow; k++) {
                    result(k,rescol) = probs[offset + k + j*nrow] * intcovar(k,i);
                }
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
    const int nrow  = d[0];
    const int ngen = d[1];
    const int npos = d[2];
    const int nint  = intcovar.cols();

    if(intcovar.rows() != nrow)
        throw std::range_error("nrow(intcovar) != nrow(probs)");

    const int ngen_result = d[1]*(nint+1); // no. cols in result
    const int recsize = nrow*ngen; // ind x geno rectangle
    const int recsize_result = nrow*ngen_result; // ind x geno rectangle in result

    NumericVector result(recsize_result*npos);

    for(int i=0; i<npos; i++) {
        // paste probs into first batch
        std::copy(probs.begin()+i*recsize,
                  probs.begin()+(i+1)*recsize,
                  result.begin()+i*recsize_result);
        for(int j=0; j<nint; j++) {
            for(int k=0; k<ngen; k++) {
                for(int s=0; s<nrow; s++)
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
    const int nrow = mat.rows();
    const int ncol = mat.cols();
    if(nrow != weights.size())
        throw std::range_error("nrow(mat) != length(weights)");

    NumericMatrix result(nrow,ncol);

    for(int j=0; j<ncol; j++)
        for(int i=0; i<nrow; i++)
            result(i,j) = mat(i,j)*weights[i];

    return result;
}

// multiply each column of a 3-dimensional array by a set of weights
// [[Rcpp::export]]
NumericVector weighted_3darray(const NumericVector& array,
                               const NumericVector& weights)
{
    const Dimension d = array.attr("dim");
    const int n = d[0];
    const int ncol = d[1]*d[2];
    if(n != weights.size())
        throw std::range_error("nrow(array) != length(weights)");

    NumericVector result(n*ncol);
    result.attr("dim") = d;

    for(int j=0, k=0; j<ncol; j++)
        for(int i=0; i<n; i++, k++)
            result[k] = array[k]*weights[i];

    return result;
}

// matrix multiplication
// [[Rcpp::export]]
NumericMatrix matrix_x_matrix(const NumericMatrix& X,
                              const NumericMatrix& Y)
{
    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd YY(as<Map<MatrixXd> >(Y));

    if(XX.cols() != YY.rows())
        throw std::range_error("ncol(X) != nrow(Y)");

    NumericMatrix result(wrap(XX * YY));
    return(result);
}

// multiply matrix by vector
// [[Rcpp::export]]
NumericVector matrix_x_vector(const NumericMatrix& X,
                              const NumericVector& y)
{
    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const VectorXd yy(as<Map<VectorXd> >(y));

    if(XX.cols() != yy.size())
        throw std::range_error("ncol(X) != length(y)");

    NumericVector result(wrap(XX * yy));
    return(result);
}

// multiply matrix by array
// [[Rcpp::export]]
NumericVector matrix_x_3darray(const NumericMatrix& X,
                               NumericVector& A)
{
    if(Rf_isNull(A.attr("dim")))
        throw std::invalid_argument("A has no dimension attribute");
    Dimension d = A.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("A should be 3-dimensional array");
    const int Xrow = X.rows();
    const int Xcol = X.cols();
    const int Arow = d[0];
    const int Acol = d[1];
    const int Apos = d[2];
    if(Xcol != Arow)
        throw std::invalid_argument("ncol(X) != nrow(A)");

    // treat as a matrix
    A.attr("dim") = Dimension(Arow, Acol*Apos);

    // cast for Eigen
    const MatrixXd XX(as<Map<MatrixXd> >(X));
    const MatrixXd AA(as<Map<MatrixXd> >(A));

    // matrix multiplication
    NumericVector result = wrap(XX*AA);
    result.attr("dim") = Dimension(Xrow, Acol, Apos);

    // fix dimension attribute of input matrix
    A.attr("dim") = d;

    return result;
}

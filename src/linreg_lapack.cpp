// linear regression via LAPACK

#include <Rcpp.h>
#include <R_ext/Lapack.h>

using namespace Rcpp;

#include "linreg_lapack.h"

// perform linear regression via LAPACK
// This does the major work; called by calc_rss_lapack or calc_resid_lapack
//
// return value contains coefficients and, for dgels, stuff that sums to RSS
// rank and used_dgelsy are needed in the functions that call this
NumericMatrix calc_regutil_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                                  int& rank, bool& used_dgelsy,
                                  const bool skip_dgels=false, const double tol=1e-10)
{
    used_dgelsy=false;
    int nrow = X.rows();
    int ncolx = X.cols(), ncoly = Y.cols();

    int minXdim = std::min(nrow, ncolx);
    int n_dwork = std::max(minXdim + std::max(minXdim, ncoly),
                           std::max(minXdim + 3*ncolx + 1, 2*minXdim*ncoly));
    std::vector<double> dwork(n_dwork);
    NumericMatrix result(nrow, ncoly);

    char notranspose='N';
    int lda=nrow, ldb=nrow, info, rank_calc;
    bool singular=false;

    std::vector<double> xx = as<std::vector<double> >(X);
    std::vector<double> yy = as<std::vector<double> >(Y);

    if(!skip_dgels) {
        F77_CALL(dgels)(&notranspose, &nrow, &ncolx, &ncoly, &xx[0],
                        &lda, &yy[0], &ldb,
                        &dwork[0], &n_dwork, &info);

        // X contains QR decomposition; if diagonal element is 0, input is rank-deficient
        rank = ncolx;
        for(int i=0, col_index=0; i<ncolx; i++, col_index += nrow) {
            if(fabs(xx[i+col_index]) < tol) {
                singular = true;
                break;
            }
        }
    }

    if(singular) {
        // restore X and Y
        std::copy(X.begin(), X.end(), xx.begin());
        std::copy(Y.begin(), Y.end(), yy.begin());
    }

    if(skip_dgels || singular) {
        std::vector<int> jpvt(ncolx);

        // use dgelsy to determine which X columns to use
        F77_CALL(dgelsy)(&nrow, &ncolx, &ncoly, &xx[0], &lda,
                         &yy[0], &ldb, &jpvt[0],
                         &tol, &rank_calc, &dwork[0], &n_dwork, &info);
        rank = rank_calc;
        used_dgelsy=true;
    }

    std::copy(yy.begin(), yy.end(), result.begin());
    return result;
}



// calculate RSS for linear regression via LAPACK
// [[Rcpp::export]]
NumericVector calc_rss_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                              const bool skip_dgels=false, const double tol=1e-10)
{
    int nrow = X.rows();
    int ncolx = X.cols(), ncoly = Y.cols();
    NumericVector rss(ncoly);
    bool used_dgelsy;
    int rank;

    // do the regression
    NumericMatrix YY = calc_regutil_lapack(X, Y, rank, used_dgelsy, skip_dgels, tol);

    if(used_dgelsy) { // X was singular so used dgelsy
        double fitted;
        for(int j=0, col_j=0; j<ncoly; j++, col_j += nrow) {
            rss[j] = 0.0;
            for(int i=0; i<nrow; i++) {

                fitted=0.0;
                for(int k=0; k<ncolx; k++) // YY contains estimated betas
                    fitted += X(i,k) * YY[k + col_j];

                double resid = Y(i,j) - fitted;
                rss[j] += resid*resid;
            }
        }
        return rss;
    }

    // calculate RSS
    for(int i=0, col_i=0; i<ncoly; i++, col_i += nrow) {
        for(int j=rank; j<nrow; j++) { // later YY's contain the residuals
            double resid = YY[j+col_i];
            rss[i] += resid*resid;
        }
    }

    return rss;
}


// calculate residuals from linear regression via LAPACK
// [[Rcpp::export]]
NumericMatrix calc_resid_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                                const bool skip_dgels=false, const double tol=1e-10)
{
    int nrow = X.rows();
    int ncolx = X.cols(), ncoly = Y.cols();
    NumericMatrix resid(nrow, ncoly);
    bool used_dgelsy;
    int rank;

    // do the regression
    NumericMatrix YY = calc_regutil_lapack(X, Y, rank, used_dgelsy, skip_dgels, tol);

    // calculate residuals
    for(int j=0, col_j=0; j<ncoly; j++, col_j += nrow) {
        for(int i=0; i<nrow; i++) {

            double fitted=0.0;
            for(int k=0; k<ncolx; k++) // YY's contain estimated betas
                fitted += X(i,k) * YY[k + col_j];

            resid(i,j) = Y(i,j) - fitted;
        }
    }

    return resid;
}

// linear regression via LAPACK

#include <Rcpp.h>
#include <R_ext/Lapack.h>

using namespace Rcpp;

#include "linreg_lapack.h"

// calculate RSS for linear regression via LAPACK
// [[Rcpp::export]]
NumericVector calc_rss_lapack(const NumericMatrix X, const NumericMatrix Y,
                              const bool skip_dgels=false, const double tol=1e-10)
{
    int nrow = X.rows();
    int ncolx = X.cols(), ncoly = Y.cols();

    int minXdim = std::min(nrow, ncolx);
    int n_dwork = std::max(minXdim + std::max(minXdim, ncoly),
                           std::max(minXdim + 3*ncolx + 1, 2*minXdim*ncoly));
    std::vector<double> dwork(n_dwork);
    NumericVector rss(ncoly);

    char notranspose='N';
    int lda=nrow, ldb=nrow, info, rank;
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
                         &tol, &rank, &dwork[0], &n_dwork, &info);


        if(rank < ncolx) {
            double fitted;
            for(int j=0, col_j=0; j<ncoly; j++, col_j += nrow) {
                rss[j] = 0.0;
                for(int i=0; i<nrow; i++) {

                    fitted=0.0;
                    for(int k=0; k<ncolx; k++) // yy's contain estimated betas
                        fitted += X(i,k) * yy[k + col_j];

                    double resid = Y(i,j) - fitted;
                    rss[j] += resid*resid;
                }
            }
            return rss;
        }
    }

    /* calculate RSS */
    for(int i=0, col_i=0; i<ncoly; i++, col_i += nrow) {
        for(int j=rank; j<nrow; j++) { // later yy's contain the residuals
            double resid = yy[j+col_i];
            rss[i] += resid*resid;
        }
    }

    return rss;
}

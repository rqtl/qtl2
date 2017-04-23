// logistic regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "binreg_eigen.h"
#include <RcppEigen.h>
#include "matrix.h"
#include "linreg_eigen.h"
#include "r_message.h"

using namespace Rcpp;
using namespace Eigen;


// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
// [[Rcpp::export]]
double calc_ll_binreg_eigenchol(const NumericMatrix& X, const NumericVector& y,
                                const int maxit=100, const double tol=1e-6)
{
    const int n_ind = y.size();
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");

    double curllik = 0;
    NumericVector pi(n_ind), wt(n_ind), nu(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind] + 0.5)/2;
        wt[ind] = sqrt(pi[ind] * (1-pi[ind]));
        nu[ind] = log(pi[ind]) - log(1-pi[ind]);
        z[ind] = nu[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
        curllik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
    }

    NumericMatrix XX = weighted_matrix(X, wt);

    bool converged=false;
    double llik;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        nu = calc_fitted_linreg_eigenchol(XX, z);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            nu[ind] /= wt[ind]; // need to divide by previous weights
            pi[ind] = exp(nu[ind])/(1.0 + exp(nu[ind]));
            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind]));
            z[ind] = nu[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
            llik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
        }

        if(fabs(llik - curllik) < tol) { // converged
            converged = true;
            break;
        }

        XX = weighted_matrix(X, wt);
        curllik = llik;

    } // end iterations

    if(!converged) r_warning("binreg didn't converge");

    return llik;
}

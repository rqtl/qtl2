// logistic regression, with weights, via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "binreg_weighted_eigen.h"
#include <RcppEigen.h>
#include "matrix.h"
#include "linreg_eigen.h"
#include "r_message.h"

using namespace Rcpp;
using namespace Eigen;


// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
// this version with weights
// [[Rcpp::export]]
double calc_ll_binreg_weighted_eigenchol(const NumericMatrix& X, const NumericVector& y,
                                         const NumericVector &weights,
                                         const int maxit=100, const double tol=1e-6)
{
    const int n_ind = y.size();
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    if(n_ind != weights.size())
        throw std::invalid_argument("nrow(X) != length(weights)");

    double curllik = 0;
    NumericVector pi(n_ind), wt(n_ind), nu(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind]*weights[ind] + 0.5)/(weights[ind] + 1.0);
        wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
        nu[ind] = log(pi[ind]) - log(1.0 - pi[ind]);
        z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
        curllik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
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
            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
            z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
            llik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
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


// logistic regression by Qr decomposition with column pivoting
// return just the log likelihood
// this version with weights
// [[Rcpp::export]]
double calc_ll_binreg_weighted_eigenqr(const NumericMatrix& X, const NumericVector& y,
                                       const NumericVector& weights,
                                       const int maxit=100, const double tol=1e-6,
                                       const double qr_tol=1e-12)
{
    const int n_ind = y.size();
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    if(n_ind != weights.size())
        throw std::invalid_argument("nrow(X) != length(weights)");

    double curllik = 0;
    NumericVector pi(n_ind), wt(n_ind), nu(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind]*weights[ind] + 0.5)/(weights[ind] + 1.0);
        wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
        nu[ind] = log(pi[ind]) - log(1.0 - pi[ind]);
        z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
        curllik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
    }

    NumericMatrix XX = weighted_matrix(X, wt);

    bool converged=false;
    double llik;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        nu = calc_fitted_linreg_eigenqr(XX, z, qr_tol);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            nu[ind] /= wt[ind]; // need to divide by previous weights
            pi[ind] = exp(nu[ind])/(1.0 + exp(nu[ind]));
            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
            z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
            llik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
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

// logistic regression
// return just the coefficients
// this version with weights
// [[Rcpp::export]]
NumericVector calc_coef_binreg_weighted_eigenqr(const NumericMatrix& X,
                                                const NumericVector& y,
                                                const NumericVector& weights,
                                                const int maxit=100,
                                                const double tol=1e-6,
                                                const double qr_tol=1e-12)
{
    List fit = fit_binreg_weighted_eigenqr(X, y, weights, false, maxit, tol, qr_tol);

    NumericVector coef = fit[2];
    return(coef);
}

// logistic regression
// return the coefficients and SEs
// this version with weights
// [[Rcpp::export]]
List calc_coefSE_binreg_weighted_eigenqr(const NumericMatrix& X,
                                         const NumericVector& y,
                                         const NumericVector & weights,
                                         const int maxit=100,
                                         const double tol=1e-6,
                                         const double qr_tol=1e-12)
{
    List fit = fit_binreg_weighted_eigenqr(X, y, weights, true, maxit, tol, qr_tol);

    NumericVector coef = fit[2];
    NumericVector SE = fit[3];

    return List::create(Named("coef") = coef,
                        Named("SE") = SE);
}

// logistic regression
// return llik, fitted probabilities, coef, SE
// this version with weights
// [[Rcpp::export]]
List fit_binreg_weighted_eigenqr(const NumericMatrix& X,
                                 const NumericVector& y,
                                 const NumericVector &weights,
                                 const bool se=true, // whether to include SEs
                                 const int maxit=100,
                                 const double tol=1e-6,
                                 const double qr_tol=1e-12)
{
    const int n_ind = y.size();
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    if(n_ind != weights.size())
        throw std::invalid_argument("nrow(X) != length(weights)");

    double curllik = 0;
    NumericVector pi(n_ind), wt(n_ind), nu(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind]*weights[ind] + 0.5)/(weights[ind] + 1.0);
        wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
        nu[ind] = log(pi[ind]) - log(1.0 - pi[ind]);
        z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
        curllik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
    }

    NumericMatrix XX = weighted_matrix(X, wt); // to store weighted matrix

    bool converged=false;
    double llik;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        nu = calc_fitted_linreg_eigenqr(XX, z, qr_tol);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            nu[ind] /= wt[ind]; // need to divide by previous weights
            pi[ind] = exp(nu[ind])/(1.0 + exp(nu[ind]));
            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind])*weights[ind]);
            z[ind] = wt[ind]*(nu[ind] + (y[ind] - pi[ind])/(pi[ind]*(1.0-pi[ind])));
            llik += (y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]))*weights[ind];
        }

        XX = weighted_matrix(X, wt);

        if(fabs(llik - curllik) < tol) { // converged
            converged = true;
            break;
        }

        curllik = llik;
    } // end iterations

    if(!converged) r_warning("binreg didn't converge");

    // now get coefficients, SEs, etc.
    List fit = fit_linreg_eigenqr(XX, z, true, qr_tol);
    NumericVector coef = fit[0];

    if(se) {
        // SE scaled by sigma; need to unscale
        NumericVector sigma = fit[4];
        NumericVector SE = fit[7];
        for(int i=0; i<SE.size(); i++) SE[i] /= sigma[0];

        return List::create(Named("log10lik") = llik,
                            Named("fitted_probs") = pi,
                            Named("coef") = coef,
                            Named("SE") = SE);
    }
    else { // no need for SEs
        return List::create(Named("log10lik") = llik,
                            Named("fitted_probs") = pi,
                            Named("coef") = coef);
    }
}

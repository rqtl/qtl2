// logistic regression via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include "binreg_eigen.h"
#include <RcppEigen.h>
#include "matrix.h"
#include "linreg_eigen.h"
#include "r_message.h" // defines RQTL2_NODEBUG + r_warning()

using namespace Rcpp;
using namespace Eigen;


// logistic regression by "LLt" Cholesky decomposition
// return just the log likelihood
//
// maxit = maximum number of iterations
// tol = tolerance for convergence
// eta_max = maximum value for linear predictor
//
// [[Rcpp::export]]
double calc_ll_binreg_eigenchol(const NumericMatrix& X, const NumericVector& y,
                                const int maxit=100, const double tol=1e-6,
                                const double eta_max=30.0)
{
    const int n_ind = y.size();
    #ifndef RQTL2_NODEBUG
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    double curllik = 0.0;
    NumericVector pi(n_ind), wt(n_ind), eta(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind] + 0.5)/2;
        wt[ind] = sqrt(pi[ind] * (1.0-pi[ind]));
        eta[ind] = log(pi[ind]) - log(1.0-pi[ind]);
        z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
        curllik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
    }

    NumericMatrix XX = weighted_matrix(X, wt);

    bool converged=false;
    double llik=0.0;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        eta = calc_fitted_linreg_eigenchol(XX, z);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            eta[ind] /= wt[ind]; // need to divide by previous weights

            // don't let eta get too large or too small
            if(eta[ind] < -eta_max) eta[ind] = -eta_max;
            else if(eta[ind] > eta_max) eta[ind] = eta_max;

            pi[ind] = exp(eta[ind])/(1.0 + exp(eta[ind]));

            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind]));
            z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
            llik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
        }

        if(fabs(llik - curllik) < tol) { // converged
            converged = true;
            break;
        }

        XX = weighted_matrix(X, wt);
        curllik = llik;

    } // end iterations

    if(!converged) r_warning("binary trait regression didn't converge: increase maxit or tol");

    return llik;
}


// logistic regression by Qr decomposition with column pivoting
// return just the log likelihood
//
// maxit = maximum number of iterations
// tol = tolerance for convergence
// qr_tol = tolerance for QR decomposition
// eta_max = maximum value for linear predictor
//
// [[Rcpp::export]]
double calc_ll_binreg_eigenqr(const NumericMatrix& X, const NumericVector& y,
                              const int maxit=100, const double tol=1e-6,
                              const double qr_tol=1e-12, const double eta_max=30.0)
{
    const int n_ind = y.size();
    #ifndef RQTL2_NODEBUG
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    double curllik = 0.0;
    NumericVector pi(n_ind), wt(n_ind), eta(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind] + 0.5)/2;
        wt[ind] = sqrt(pi[ind] * (1.0-pi[ind]));
        eta[ind] = log(pi[ind]) - log(1.0-pi[ind]);
        z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
        curllik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
    }

    NumericMatrix XX = weighted_matrix(X, wt);

    bool converged=false;
    double llik=0.0;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        eta = calc_fitted_linreg_eigenqr(XX, z, qr_tol);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            eta[ind] /= wt[ind]; // need to divide by previous weights

            // don't let eta get too large or too small
            if(eta[ind] < -eta_max) eta[ind] = -eta_max;
            else if(eta[ind] > eta_max) eta[ind] = eta_max;

            pi[ind] = exp(eta[ind])/(1.0 + exp(eta[ind]));

            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind]));
            z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];
            llik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
        }

        if(fabs(llik - curllik) < tol) { // converged
            converged = true;
            break;
        }

        XX = weighted_matrix(X, wt);
        curllik = llik;

    } // end iterations

    if(!converged) r_warning("binary trait regression didn't converge: increase maxit or tol");

    return llik;
}

// logistic regression
// return just the coefficients
// [[Rcpp::export]]
NumericVector calc_coef_binreg_eigenqr(const NumericMatrix& X,
                                       const NumericVector& y,
                                       const int maxit=100,       // max iterations
                                       const double tol=1e-6,     // tolerance for convergence
                                       const double qr_tol=1e-12, // tolerance for QR decomp
                                       const double eta_max=30.0)  // max value for linear predictor
{
    List fit = fit_binreg_eigenqr(X, y, false, maxit, tol, qr_tol, eta_max);

    NumericVector coef = fit[2];

    return coef;
}

// logistic regression
// return the coefficients and SEs
// [[Rcpp::export]]
List calc_coefSE_binreg_eigenqr(const NumericMatrix& X,
                                const NumericVector& y,
                                const int maxit=100,       // max iterations
                                const double tol=1e-6,     // tolerance for convergence
                                const double qr_tol=1e-12, // tolerance for QR decomp
                                const double eta_max=30.0)  // max value for linear predictor
{
    List fit = fit_binreg_eigenqr(X, y, true, maxit, tol, qr_tol, eta_max);

    NumericVector coef = fit[2];
    NumericVector SE = fit[3];

    return List::create(Named("coef") = coef,
                        Named("SE") = SE);
}

// logistic regression
// return llik, fitted probabilities, coef, SE
// [[Rcpp::export]]
List fit_binreg_eigenqr(const NumericMatrix& X,
                        const NumericVector& y,
                        const bool se=true,        // whether to include SEs
                        const int maxit=100,       // max iterations
                        const double tol=1e-6,     // tolerance for convergence
                        const double qr_tol=1e-12, // tolerance for QR decomp
                        const double eta_max=30.0)  // max value for linear predictor
{
    const int n_ind = y.size();
    #ifndef RQTL2_NODEBUG
    if(n_ind != X.rows())
        throw std::invalid_argument("nrow(X) != length(y)");
    #endif

    double curllik = 0.0;
    NumericVector pi(n_ind), wt(n_ind), eta(n_ind), z(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        pi[ind] = (y[ind] + 0.5)/2.0;
        wt[ind] = sqrt(pi[ind] * (1.0-pi[ind]));
        eta[ind] = log(pi[ind]) - log(1.0-pi[ind]);
        z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];

        curllik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
    }

    NumericMatrix XX = weighted_matrix(X, wt); // to store weighted matrix

    bool converged=false;
    double llik=0.0;

    for(int it=0; it<maxit; it++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // fitted values using weighted XX; will need to divide by previous weights
        eta = calc_fitted_linreg_eigenqr(XX, z, qr_tol);

        llik = 0.0;
        for(int ind=0; ind<n_ind; ind++) {
            eta[ind] /= wt[ind]; // need to divide by previous weights

            // don't let eta get too large or too small
            if(eta[ind] < -eta_max) eta[ind] = -eta_max;
            else if(eta[ind] > eta_max) eta[ind] = eta_max;

            pi[ind] = exp(eta[ind])/(1.0 + exp(eta[ind]));

            wt[ind] = sqrt(pi[ind] * (1.0 - pi[ind]));
            z[ind] = eta[ind]*wt[ind] + (y[ind] - pi[ind])/wt[ind];

            llik += y[ind] * log10(pi[ind]) + (1.0-y[ind])*log10(1.0-pi[ind]);
        }

        XX = weighted_matrix(X, wt);

        if(fabs(llik - curllik) < tol) { // converged
            converged = true;
            break;
        }

        curllik = llik;
    } // end iterations

    if(!converged) r_warning("binary trait regression didn't converge: increase maxit or tol");

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

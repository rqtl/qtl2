// HMM utility functions

#include "hmm_util.h"
#include <math.h>
#include <Rcpp.h>

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
// [[Rcpp::export]]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    // if both -Inf, return -Inf
    if(Rcpp::traits::is_infinite<REALSXP>(a) &&
       Rcpp::traits::is_infinite<REALSXP>(b) && a < 0 && b < 0) return a;

    if(b > a + tol) return b;
    else if(a > b + tol) return a;
    else return a + log1p(exp(b-a));
}

// Calculate subtractlog(a,b) = log[exp(a) - exp(b)]
// [[Rcpp::export]]
double subtractlog(const double a, const double b)
{
    const double tol=200.0;

    if(a > b + tol) return a;
    else return a + log1p(-exp(b-a));
}

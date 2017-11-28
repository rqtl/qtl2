// HMM utility functions

#include "hmm_util.h"
#include <math.h>
#include <Rcpp.h>

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
// [[Rcpp::export]]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    if(b > a + tol) return b;
    else if(a > b + tol) return a;
    else return a + log1p(exp(b-a));
}

// Calculate subtrlog(a,b) = log[exp(a) - exp(b)]
// [[Rcpp::export]]
double subtrlog(const double a, const double b)
{
    const double tol=200.0;

    if(a > b + tol) return a;
    else return a + log1p(-exp(b-a));
}

// main HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm.h"

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
// [[Rcpp::export]]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    if(b > a + tol) return b;
    else if(a > b + tol) return a;
    else return a + log1p(exp(b-a));
}

// Utility functions

<<<<<<< HEAD:src/hmm_util.cpp
#include "hmm_util.h"
#include <math.h>
#include <Rcpp.h>
=======
#include "util.h"
#include <math.h>
>>>>>>> qtl2scan/master:src/util.cpp

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    if(b > a + tol) return b;
    else if(a > b + tol) return a;
    else return a + log1p(exp(b-a));
}

<<<<<<< HEAD:src/hmm_util.cpp
// Calculate  subtrlog(a,b) = log[exp(a) - exp(b)]
// [[Rcpp::export]]
double subtrlog(const double a, const double b)
=======
// Calculate subtractlog(a,b) = log[exp(a) - exp(b)]
double subtractlog(const double a, const double b)
>>>>>>> qtl2scan/master:src/util.cpp
{
    const double tol=200.0;

    if(a > b + tol) return a;
    else return a + log1p(-exp(b-a));
}

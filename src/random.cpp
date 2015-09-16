// random number generation

#include <Rcpp.h>
using namespace Rcpp;

#include "random.h"

// sample random integer from 0, 1, 2, ..., n-1 with probability p[0], p[1], ...
int sample_int(NumericVector probs)
{
    int n=probs.size();
    int result;

    double u = R::runif(0.0, 1.0);

    for(result=0; result < n; result++) {
        if(u <= probs[result]) return result;
        u -= probs[result];
    }

    return NA_INTEGER;
}

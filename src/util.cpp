// Utility functions

#include "util.h"
#include <math.h>

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    if(b > a + tol) return b;
    else if(a > b + tol) return a;
    else return a + log1p(exp(b-a));
}

// Calculate subtractlog(a,b) = log[exp(a) - exp(b)]
double subtractlog(const double a, const double b)
{
    const double tol=200.0;

    if(a > b + tol) return a;
    else return a + log1p(-exp(b-a));
}

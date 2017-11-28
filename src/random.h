// random number generation
#ifndef RANDOM_H
#define RANDOM_H

#include <Rcpp.h>

// sample random integer from 0, 1, 2, ..., n-1 with probability p[0], p[1], ...
int sample_int(Rcpp::NumericVector probs);

// sample random integer from 0, 1, 2, ..., n-1, with equal probabilities
int sample_int(int n);

#endif // RANDOM_H

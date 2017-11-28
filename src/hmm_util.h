// HMM utility functions
#ifndef HMM_UTIL_H
#define HMM_UTIL_H

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
double addlog(const double a, const double b);

// Calculate  subtrlog(a,b) = log[exp(a) - exp(b)]
double subtrlog(const double a, const double b);

#endif // HMM_UTIL_H

// 1-d optimization by Brent's method
#ifndef BRENT_FMIN_H
#define BRENT_FMIN_H

/***********************************************************************
 * Taken from R ver 3.2.2
 * R-3.2.2/src/library/stats/src/optimize.c
 ***********************************************************************/

// see brent_fmin.cpp

double lmm_Brent_fmin(double ax, double bx,
                      double f(const double x, struct calcLL_args *info),
                      calcLL_args *info, double tol);

#endif // BRENT_FMIN_H

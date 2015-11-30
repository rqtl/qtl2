// generic linear regression functions
#ifndef LINREG_H
#define LINREG_H

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                              const double tol);

// Calculate matrix of residuals from linear regression of Y on X
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y,
                                const double tol);

// use calc_resid_linreg for a 3-dim array
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P,
                                   const double tol);

#endif // LINREG_H

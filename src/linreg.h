// generic linear regression functions

// Calculate vector of residual sum of squares (RSS) from linear regression of Y vs X
NumericVector calc_rss_linreg(const NumericMatrix& X, const NumericMatrix& Y);

// Calculate matrix of residuals from linear regression of Y on X
NumericMatrix calc_resid_linreg(const NumericMatrix& X, const NumericMatrix& Y);

// use calc_resid_linreg for a 3-dim array
NumericVector calc_resid_linreg_3d(const NumericMatrix& X, const NumericVector& P);

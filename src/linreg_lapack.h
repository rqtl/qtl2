// calculate RSS for linear regression via LAPACK
NumericVector calc_rss_lapack(const NumericMatrix X, const NumericMatrix Y,
                              bool skip_dgels, double tol);

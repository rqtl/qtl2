// calculate RSS for linear regression via LAPACK
NumericVector calc_rss_lapack(const NumericMatrix X, const NumericMatrix Y,
                              const bool skip_dgels, const double tol);

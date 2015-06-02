// linear regression via LAPACK

// perform linear regression via LAPACK
// This does the major work; called by calc_rss_lapack or calc_resid_lapack
//
// return value contains coefficients and, for dgels, stuff that sums to RSS
// rank and used_dgelsy are needed in the functions that call this
NumericMatrix calc_regutil_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                                  int& rank, bool& used_dgelsy,
                                  const bool skip_dgels, const double tol);

// calculate RSS for linear regression via LAPACK
NumericVector calc_rss_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                              const bool skip_dgels, const double tol);

// calculate residuals from linear regression via LAPACK
NumericMatrix calc_resid_lapack(const NumericMatrix& X, const NumericMatrix& Y,
                                const bool skip_dgels, const double tol);

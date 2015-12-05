// get X chromosome covariates
#ifndef GET_X_COVAR_H
#define GET_X_COVAR_H

#include <Rcpp.h>

Rcpp::NumericMatrix get_x_covar(const Rcpp::String& crosstype,
                                const Rcpp::LogicalVector& is_female, // length n_ind
                                const Rcpp::IntegerMatrix& cross_info); // columns are individuals

#endif // GET_X_COVAR_H

// get X chromosome covariates
#ifndef GET_X_COVAR_H
#define GET_X_COVAR_H

NumericMatrix get_x_covar(const String& crosstype,
                          const LogicalVector& is_female, // length n_ind
                          const IntegerMatrix& cross_info); // columns are individuals

#endif // GET_X_COVAR_H

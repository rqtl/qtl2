// get X chromosome covariates

#include <Rcpp.h>
#include "cross.h"
#include "get_x_covar.h"

// [[Rcpp::export(".get_x_covar")]]
NumericMatrix get_x_covar(const String& crosstype,
                          const LogicalVector& is_female, // length n_ind
                          const IntegerMatrix& cross_info) // columns are individuals
{
    QTLCross* cross = QTLCross::Create(crosstype);

    return cross->get_x_covar(is_female, cross_info);
}

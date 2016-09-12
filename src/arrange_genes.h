// arrange_genes, helper function for plot_genes
#ifndef QTL2PLOT_ARRANGE_GENES
#define QTL2PLOT_ARRANGE_GENES

// take vectors (start, stop) and return y values in {1, 2, ...}
// so that the genes won't overlap

#include <Rcpp.h>

Rcpp::IntegerVector arrange_genes(const Rcpp::NumericVector& start,
                                  const Rcpp::NumericVector& end);

#endif // QTL2PLOT_ARRANGE_GENES

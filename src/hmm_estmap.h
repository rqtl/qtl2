// re-estimate inter-marker recombination fractions
#ifndef HMM_ESTMAP_H
#define HMM_ESTMAP_H

#include <Rcpp.h>

Rcpp::NumericVector est_map(const Rcpp::String& crosstype,
                            const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                            const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                            const bool is_X_chr,
                            const Rcpp::LogicalVector& is_female,
                            const Rcpp::IntegerMatrix& cross_info,
                            const Rcpp::NumericVector& rec_frac,
                            const double error_prob,
                            const int max_iterations,
                            const double tol,
                            const bool verbose);

#endif // HMM_ESTMAP_H

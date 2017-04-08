// re-estimate inter-marker recombination fractions
// this set of versions is tailored somewhat to cross type, to be faster
// est_map2_simple: same transition matrix for all individuals
//                  (bc-type, f2, riself4)
// est_map2_forder: need to deal with founder order in transition matrix
//                  (riself8, riself16)
// est_map2_ngen:   need separate transition matrices for each unique value of
//                  number of generations (ail, do) [not yet implemented]
#ifndef HMM_ESTMAP2_H
#define HMM_ESTMAP2_H

#include <Rcpp.h>

// re-estimate inter-marker recombination fractions
Rcpp::NumericVector est_map2(const Rcpp::String& crosstype,
                             const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                             const Rcpp::IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                             const bool is_X_chr,
                             const Rcpp::LogicalVector& is_female,
                             const Rcpp::IntegerMatrix& cross_info,
                             const Rcpp::IntegerVector& cross_group,
                             const Rcpp::IntegerVector& unique_cross_group,
                             const Rcpp::NumericVector& rec_frac,
                             const double error_prob,
                             const int max_iterations,
                             const double tol,
                             const bool verbose);


// actually just using the low-mem approach here for now
Rcpp::NumericVector est_map2_simple(const Rcpp::String crosstype,
                                    const Rcpp::IntegerMatrix& genotypes,
                                    const Rcpp::IntegerMatrix& founder_geno,
                                    const bool is_X_chr,
                                    const Rcpp::LogicalVector& is_female,
                                    const Rcpp::IntegerMatrix& cross_info,
                                    const Rcpp::IntegerVector& cross_group,
                                    const Rcpp::IntegerVector& unique_cross_group,
                                    const Rcpp::NumericVector& rec_frac,
                                    const double error_prob,
                                    const int max_iterations,
                                    const double tol,
                                    const bool verbose);

// Need same set of possible genotypes for all individuals,
// and same basic structure for transition matrix, but reorder transition matrix by founder order
// (for riself8 and riself16)
Rcpp::NumericVector est_map2_founderorder(const Rcpp::String crosstype,
                                          const Rcpp::IntegerMatrix& genotypes,
                                          const Rcpp::IntegerMatrix& founder_geno,
                                          const bool is_X_chr,
                                          const Rcpp::LogicalVector& is_female,
                                          const Rcpp::IntegerMatrix& cross_info,
                                          const Rcpp::IntegerVector& cross_group,
                                          const Rcpp::IntegerVector& unique_cross_group,
                                          const Rcpp::NumericVector& rec_frac,
                                          const double error_prob,
                                          const int max_iterations,
                                          const double tol,
                                          const bool verbose);

#endif // HMM_ESTMAP2_H

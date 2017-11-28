// functions to check QTL cross data/information
#ifndef CHECK_CROSS_H
#define CHECK_CROSS_H

#include <Rcpp.h>

// check if a cross type is supported
bool crosstype_supported(const Rcpp::String& crosstype);

// count inconsistencies in marker data
Rcpp::IntegerVector count_invalid_genotypes(const Rcpp::String& crosstype,
                                            const Rcpp::IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                            const bool& is_X_chr,
                                            const Rcpp::LogicalVector& is_female,
                                            const Rcpp::IntegerMatrix& cross_info); // columns are individuals

// check cross info
bool check_crossinfo(const Rcpp::String& crosstype,
                     const Rcpp::IntegerMatrix& cross_info,
                     const bool any_x_chr);

// check sex
bool check_is_female_vector(const Rcpp::String& crosstype,
                            const Rcpp::LogicalVector& is_female,
                            const bool any_x_chr);

// check if X chr can be handled
bool check_handle_x_chr(const Rcpp::String& crosstype,
                        const bool any_x_chr);

#endif // CHECK_CROSS_H

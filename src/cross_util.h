// cross utility functions
#ifndef CROSS_UTIL_H
#define CROSS_UTIL_H

#include <Rcpp.h>

// allele pair -> genotype code (for multi-parent crosses with heterozygosity)
int mpp_encode_alleles(const int allele1, const int allele2,
                       const int n_alleles, const bool phase_known);

// genotype code -> allele pair (for multi-parent crosses with heterozygosity)
Rcpp::IntegerVector mpp_decode_geno(const int true_gen, const int n_alleles,
                                    const bool phase_known);

// is heterozygous? (for multi-parent crosses with heterozygosity)
bool mpp_is_het(const int true_gen, const int n_alleles, const bool phase_known);

// geno_names from allele names
const std::vector<std::string> mpp_geno_names(const std::vector<std::string> alleles,
                                              const bool is_x_chr);

#endif // CROSS_UTIL_H

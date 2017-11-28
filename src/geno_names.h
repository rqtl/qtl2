// genotype names from alleles
#ifndef GENO_NAMES_H
#define GENO_NAMES_H

#include <Rcpp.h>

std::vector<std::string> geno_names(const Rcpp::String& crosstype,
                                    const std::vector<std::string> alleles,
                                    const bool is_x_chr);

#endif // GENO_NAMES_H

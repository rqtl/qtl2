// genotype names from alleles

#include "geno_names.h"
#include <Rcpp.h>
#include "cross.h"

// [[Rcpp::export]]
std::vector<std::string> geno_names(const String& crosstype,
                                    const std::vector<std::string> alleles,
                                    const bool is_x_chr)
{
    QTLCross* cross = QTLCross::Create(crosstype);
    std::vector<std::string> result = cross->geno_names(alleles, is_x_chr);

    delete cross;

    return result;
}

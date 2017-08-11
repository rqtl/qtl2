// haploid QTLCross class (for HMM)

#include "cross_haploid.h"
#include <Rcpp.h>
#include "cross.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

// geno_names from allele names
const std::vector<std::string> HAPLOID::geno_names(const std::vector<std::string> alleles,
                                                   const bool is_x_chr)
{
    if(alleles.size() < 2)
        throw std::range_error("alleles must have length 2");

    std::vector<std::string> result(2);
    result[0] = alleles[0];
    result[1] = alleles[1];
    return result;
}

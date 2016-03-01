// Converting genotype/allele probabilities to SNP probs

#include "snpprobs.h"
#include <exception>
#include <Rcpp.h>
using namespace Rcpp;


// calculate strain distribution pattern (SDP) from
// SNP genotypes for a set of strains
//
// Input is a marker x strain matrix of genotypes
// 0 = homozygous AA, 1 = homozygous BB
//
// [[Rcpp::export]]
IntegerVector calc_sdp(const IntegerMatrix& geno)
{
    const unsigned int n_mar = geno.rows();
    const unsigned int n_str = geno.cols();
    if(n_str < 2)
        throw std::invalid_argument("Need genotypes on >= 2 strains");

    IntegerVector result(n_mar);
    for(unsigned int i=0; i<n_mar; i++) {
        for(unsigned int j=0; j<n_str; j++) {
            result[i] += geno(i,j)*(1 << j);
        }
    }

    return result;
}

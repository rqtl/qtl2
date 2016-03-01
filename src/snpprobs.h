// Converting genotype/allele probabilities to SNP probs
#ifndef SNPPROBS_H
#define SNPPROBS_H

#include "snpprobs.h"
#include <exception>
#include <Rcpp.h>
using namespace Rcpp;


// calculate strain distribution pattern (SDP) from
// SNP genotypes for a set of strains
//
// Input is a marker x strain matrix of genotypes
// 0 = homozygous AA, 1 = homozygous BB
Rcpp::IntegerVector calc_sdp(Rcpp::IntegerMatrix geno);

#endif // SNPPROBS_H

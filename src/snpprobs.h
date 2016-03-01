// Converting genotype/allele probabilities to SNP probs
#ifndef SNPPROBS_H
#define SNPPROBS_H

#include <Rcpp.h>

// calculate strain distribution pattern (SDP) from
// SNP genotypes for a set of strains
//
// Input is a marker x strain matrix of genotypes
// 0 = homozygous AA, 1 = homozygous BB
Rcpp::IntegerVector calc_sdp(const Rcpp::IntegerMatrix& geno);

// convert allele probabilities into SNP probabilities
//
// alleleprob = individual x allele x position
// map = locations of alleleprob positions
// sdp = vector of strain distribution patterns
// interval = map interval containing snp
// on_map = logical vector indicating snp is at left endpoint of interval
Rcpp::NumericVector alleleprob_to_snpprob(Rcpp::NumericVector alleleprob,
                                          Rcpp::IntegerVector sdp,
                                          Rcpp::IntegerVector interval,
                                          Rcpp::LogicalVector on_map);


// convert genotype columns to SNP columns
//
// n_str     Number of strains
//    (so n_str*(n_str+1)/2 columns)
// sdp       Strain distribution pattern for SNP
Rcpp::IntegerVector genocol_to_snpcol(const int n_str, const int sdp);

// convert genotype probabilities into SNP probabilities
//
// genoprob = individual x genotype x position
// sdp = vector of strain distribution patterns
// interval = map interval containing snp
// on_map = logical vector indicating snp is at left endpoint of interval
Rcpp::NumericVector genoprob_to_snpprob(Rcpp::NumericVector genoprob,
                                        Rcpp::IntegerVector sdp,
                                        Rcpp::IntegerVector interval,
                                        Rcpp::LogicalVector on_map);


// convert X genotype columns to SNP columns
//
// n_str     Number of strains
//    (so n_str + n_str*(n_str+1)/2 columns)
// sdp       Strain distribution pattern for SNP
Rcpp::IntegerVector genocol_to_snpcol(const int n_str, const int sdp);

// convert X chr genotype probabilities into SNP probabilities
//
// here the genotypes are the 36 female genotypes followed by the 8 male genotypes
//
// genoprob = individual x genotype x position
// sdp = vector of strain distribution patterns
// interval = map interval containing snp
// on_map = logical vector indicating snp is at left endpoint of interval
Rcpp::NumericVector Xgenoprob_to_snpprob(Rcpp::NumericVector genoprob,
                                         Rcpp::IntegerVector sdp,
                                         Rcpp::IntegerVector interval,
                                         Rcpp::LogicalVector on_map);

#endif // SNPPROBS_H

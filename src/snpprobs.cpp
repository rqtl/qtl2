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

// convert allele probabilities into SNP probabilities
//
// alleleprob = individual x allele x position
// sdp = vector of strain distribution patterns
// interval = map interval containing snp
// on_map = logical vector indicating snp is at left endpoint of interval
//
// [[Rcpp::export]]
NumericVector alleleprob_to_snpprob(NumericVector alleleprob,
                                    IntegerVector sdp,
                                    IntegerVector interval,
                                    LogicalVector on_map)
{
    const IntegerVector& d = alleleprob.attr("dim");
    const unsigned int n_ind = d[0];
    const unsigned int n_str = d[1];
    const unsigned int n_pos = d[2];
    const unsigned int n_snp = sdp.size();
    if(n_snp != interval.size())
        throw std::invalid_argument("length(sdp) != length(interval)");
    if(n_snp != on_map.size())
        throw std::invalid_argument("length(sdp) != length(on_map)");

    NumericVector result(n_ind*2*n_snp);
    result.attr("dim") = Dimension(n_ind, 2, n_snp);

    // check that the interval and SDP values are okay
    for(unsigned int i=0; i<n_snp; i++) {
        if(interval[i] < 0 || interval[i] > n_pos-1 ||
           (interval[i] == n_pos-1 && !on_map[i]))
            throw std::invalid_argument("snp outside of map range");
        if(sdp[i] < 1 || sdp[i] > (1 << n_str)-1)
            throw std::invalid_argument("SDP out of range");
    }

    for(unsigned int snp=0; snp<n_snp; snp++) {
        for(unsigned int strain=0; strain < n_str; strain++) {

            unsigned int allele = ((sdp[snp] & (1 << strain)) != 0); // 0 or 1 SNP genotype for this strain

            unsigned int result_offset = allele*n_ind + (snp*n_ind*2);
            unsigned int input_offset = strain*n_ind + (interval[snp]*n_ind*n_str);
            unsigned int next_on_map = input_offset + n_ind*n_str;
            for(unsigned int ind=0; ind<n_ind; ind++) {
                if(on_map[snp]) {
                    result[ind + result_offset] += alleleprob[ind + input_offset];
                }
                else {
                    result[ind + result_offset] += (alleleprob[ind + input_offset] +
                                                    alleleprob[ind + next_on_map])/2.0;
                }
            } // loop over individuals
        }
    } // loop over snps

    return result;
}

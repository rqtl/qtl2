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
// [[Rcpp::export(".alleleprob_to_snpprob")]]
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
    if(n_str < 3) // not meaningful for <3 strains
        throw std::invalid_argument("meaningful only with >= 3 strains");

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

// convert genotype columns to SNP columns
//
// n_str     Number of strains
//    (so n_str*(n_str+1)/2 columns)
// sdp       Strain distribution pattern for SNP
//
// [[Rcpp::export]]
IntegerVector genocol_to_snpcol(const int n_str, const int sdp)
{
    const unsigned int n_gen = n_str*(n_str+1)/2;
    if(sdp < 1 || sdp > (1 << n_str)-1)
        throw std::invalid_argument("SDP out of range");

    IntegerVector result(n_gen);

    for(int a1=0, g=0; a1<n_str; a1++) {
        for(int a2=0; a2<=a1; a2++, g++) {
            // a1,a2 are the alleles
            // g is the corresponding genotype code

            // snp alleles for these two alleles
            int snp1 = ((sdp & (1 << a1)) != 0);
            int snp2 = ((sdp & (1 << a2)) != 0);

            if(snp1==0 && snp2==0) { // hom 00
                result[g] = 0;
            }
            else if(snp1==1 && snp2==1) { // hom 11
                result[g] = 2;
            }
            else { // het
                result[g] = 1;
            }
        }
    }

    return result;
}

// convert genotype probabilities into SNP probabilities
//
// genoprob = individual x genotype x position
// sdp = vector of strain distribution patterns
// interval = map interval containing snp
// on_map = logical vector indicating snp is at left endpoint of interval
//
// [[Rcpp::export(".genoprob_to_snpprob")]]
NumericVector genoprob_to_snpprob(NumericVector genoprob,
                                  IntegerVector sdp,
                                  IntegerVector interval,
                                  LogicalVector on_map)
{
    const IntegerVector& d = genoprob.attr("dim");
    const unsigned int n_ind = d[0];
    const unsigned int n_gen = d[1];
    const unsigned int n_str = (sqrt(8*n_gen + 1) - 1)/2;
    if(n_gen != n_str*(n_str+1)/2)
        throw std::invalid_argument("n_gen must == n(n+1)/2 for some n");
    const unsigned int n_pos = d[2];
    const unsigned int n_snp = sdp.size();
    if(n_snp != interval.size())
        throw std::invalid_argument("length(sdp) != length(interval)");
    if(n_snp != on_map.size())
        throw std::invalid_argument("length(sdp) != length(on_map)");
    if(n_str < 3) // not meaningful for <3 strains
        throw std::invalid_argument("meaningful only with >= 3 strains");

    NumericVector result(n_ind*3*n_snp); // 3 is the number of SNP genotypes (AA,AB,BB)
    result.attr("dim") = Dimension(n_ind, 3, n_snp);

    // check that the interval and SDP values are okay
    for(unsigned int i=0; i<n_snp; i++) {
        if(interval[i] < 0 || interval[i] > n_pos-1 ||
           (interval[i] == n_pos-1 && !on_map[i]))
            throw std::invalid_argument("snp outside of map range");
        if(sdp[i] < 1 || sdp[i] > (1 << n_str)-1)
            throw std::invalid_argument("SDP out of range");
    }

    for(unsigned int snp=0; snp<n_snp; snp++) {
        IntegerVector snpcol = genocol_to_snpcol(n_str, sdp[snp]);

        for(unsigned int g=0; g<n_gen; g++) {
            unsigned int result_offset = snpcol[g]*n_ind + (snp*n_ind*3);
            unsigned int input_offset = g*n_ind + (interval[snp]*n_ind*n_gen);
            unsigned int next_on_map = input_offset + n_ind*n_gen;
            for(unsigned int ind=0; ind<n_ind; ind++) {
                if(on_map[snp]) {
                    result[ind + result_offset] += genoprob[ind + input_offset];
                }
                else {
                    result[ind + result_offset] += (genoprob[ind + input_offset] +
                                                    genoprob[ind + next_on_map])/2.0;
                }
            } // loop over individuals
        }
    } // loop over snps

    return result;
}

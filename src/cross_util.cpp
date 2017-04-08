// cross utility functions

#include "cross_util.h"
#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// allele pair -> genotype code (for multi-parent crosses with heterozygosity)
// [[Rcpp::export]]
int mpp_encode_alleles(const int allele1, const int allele2,
                       const int n_alleles, const bool phase_known)
{
    if(phase_known) {
        const int m = std::max(allele1, allele2);
        const int d = abs(allele1 - allele2);

        if(allele1 <= allele2)
            return (int)round(R::choose((double)(m+1), 2.0) - d);
        else
            return (int)round(R::choose((double)(m), 2.0) - d + 1 +
                              R::choose((double)(n_alleles+1), 2.0));
    }
    else {
        const int m = std::max(allele1, allele2);
        const int d = abs(allele1 - allele2);
        return (int)round(R::choose((double)(m+1), 2.0) - d);
    }
}


// genotype code -> allele pair (for multi-parent crosses with heterozygosity)
// [[Rcpp::export]]
IntegerVector mpp_decode_geno(const int true_gen,
                              const int n_alleles, const bool phase_known)
{
    IntegerVector result(2);

    if(phase_known) {
        // number of phase-known genotypes
        const int n_puk_geno = n_alleles + (int)round(R::choose((double)n_alleles, 2.0));

        #ifndef NDEBUG
        const int n_geno = n_alleles * n_alleles;
        if(true_gen < 0 || true_gen > n_geno)
            throw std::range_error("genotype value not allowed");
        #endif

        if(true_gen <= n_puk_geno) {
            int last_max = 0;
            for(int i=1; i<=n_alleles; i++) {
                if(true_gen <= last_max+i) {
                    result[1] = i;
                    result[0] = true_gen - last_max;
                    return result;
                }
                last_max += i;
            }
        }
        else {
            int g = true_gen - n_puk_geno;
            int last_max = 0;
            for(int i=1; i<=n_alleles-1; i++) {
                if(g <= last_max+i) {
                    result[0] = i+1;
                    result[1] = g-last_max;
                    return result;
                }
                last_max += i;
            }
        }

        result[0] = NA_INTEGER;
        result[1] = NA_INTEGER;
        return result;
    }
    else {
        // number of phase-unknown genotypes
        #ifndef NDEBUG
        const int n_geno = n_alleles + (int)round(R::choose((double)n_alleles, 2.0));
        if(true_gen < 0 || true_gen > n_geno)
            throw std::range_error("genotype value not allowed");
        #endif

        int last_max = 0;
        for(int i=1; i<=n_alleles; i++) {
            if(true_gen <= last_max+i) {
                result[1] = i;
                result[0] = true_gen - last_max;
                return(result);
            }
            last_max += i;
        }

        result[0] = NA_INTEGER;
        result[1] = NA_INTEGER;
        return result;
    }
}

// is heterozygous? (for multi-parent crosses with heterozygosity)
// [[Rcpp::export]]
bool mpp_is_het(const int true_gen, const int n_alleles,
                const bool phase_known)
{
    IntegerVector alleles = mpp_decode_geno(true_gen, n_alleles, phase_known);
    if(alleles[0] == alleles[1]) return(false);
    return(true);
}


// geno_names from allele names
// [[Rcpp::export]]
const std::vector<std::string> mpp_geno_names(const std::vector<std::string> alleles,
                                              const bool is_x_chr)
{
    const int n_alleles = alleles.size();
    const int n_geno = n_alleles + (int)round(R::choose((double)n_alleles, 2.0));

    if(is_x_chr) {
        std::vector<std::string> result(n_geno + n_alleles);
        for(int i=0; i<n_geno; i++) {
            IntegerVector allele_int = mpp_decode_geno(i+1, n_alleles, false);
            result[i] = alleles[allele_int[0]-1] + alleles[allele_int[1]-1];
        }
        for(int i=0; i<n_alleles; i++) {
            result[n_geno+i] = alleles[i] + "Y";
        }

        return result;
    }
    else {
        std::vector<std::string> result(n_geno);
        for(int i=0; i<n_geno; i++) {
            IntegerVector allele_int = mpp_decode_geno(i+1, n_alleles, false);
            result[i] = alleles[allele_int[0]-1] + alleles[allele_int[1]-1];
        }
        return result;
    }
}

// For vector of founder indices (cross info), create the reverse index
// with the inverted indices starting at 0
// (2,3,1,4) -> (2,0,1,3)
// (2,4,3,1) -> (3,0,2,1)
// (7,8,3,5,4,1,6,2) -> (5,7,2,4,3,6,0,1)
//
// [[Rcpp::export]]
IntegerVector invert_founder_index(IntegerVector cross_info)
{
    const int n = cross_info.size();
    IntegerVector result(n);

    for(int i=0; i<n; i++) {
        const int f = cross_info[i];
        #ifndef NDEBUG
        if(f < 1 || f > n)
            throw std::range_error("cross_info has values out of range");
        #endif
        result[f-1] = i;
    }

    return(result);
}

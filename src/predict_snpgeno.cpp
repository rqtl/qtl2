// get predicted SNP genotypes from inferred genotypes + founder genotypes

#include "predict_snpgeno.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".predict_snpgeno")]]
IntegerMatrix predict_snpgeno(const IntegerMatrix& allele1,
                              const IntegerMatrix& allele2,
                              const IntegerMatrix& founder_geno)
{
    const int n_ind = allele1.rows();
    const int n_mar = allele1.cols();
    const int n_founders = founder_geno.rows();
    if(n_ind != allele2.rows())
        throw std::invalid_argument("nrow(allele1) != nrow(allele2)");
    if(n_mar != allele2.cols())
        throw std::invalid_argument("ncol(allele1) != ncol(allele2)");
    if(n_mar != founder_geno.cols())
        throw std::invalid_argument("ncol(allele1) != ncol(founder_geno)");

    IntegerMatrix result(n_ind, n_mar);
    for(int ind=0; ind<n_ind; ind++) {
        for(int mar=0; mar<n_mar; mar++) {
            if(!IntegerVector::is_na(allele1(ind,mar)) &&
               !IntegerVector::is_na(allele2(ind,mar)) &&
               founder_geno(allele1(ind,mar)-1,mar)!=0 &&
               founder_geno(allele2(ind,mar)-1,mar)!=0 &&
               allele1(ind,mar) <= n_founders &&
               allele2(ind,mar) <= n_founders) {

                result(ind,mar) = (founder_geno(allele1(ind,mar)-1,mar) - 1)/2 +
                    (founder_geno(allele2(ind,mar)-1,mar)-1)/2 + 1;

            } else {
                result(ind,mar) = NA_INTEGER;
            }
        }
    }

    return result;
}

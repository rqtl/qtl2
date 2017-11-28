// calculate matrix of counts of genotype matches for pairs of individuals

#include "compare_geno.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(".compare_geno")]]
IntegerMatrix compare_geno(const IntegerMatrix& geno) // matrix n_mar x n_ind (transposed of normal)
{
    const int n_mar = geno.rows();
    const int n_ind = geno.cols();

    IntegerMatrix result(n_ind, n_ind);

    for(int ind_i=0; ind_i<n_ind; ++ind_i) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        // diagonal is number of genotypes for individual i
        int n_typed=0;
        for(int mar=0; mar<n_mar; mar++)
            if(geno(mar, ind_i) > 0) n_typed++;
        result(ind_i, ind_i) = n_typed;

        for(int ind_j=ind_i+1; ind_j<n_ind; ind_j++) {

            int n_typed = 0.0;
            int n_matches = 0.0;
            for(int mar=0; mar<n_mar; mar++) {
                if(geno(mar, ind_i) > 0 && geno(mar, ind_j)>0) {
                    n_typed++;
                    if(geno(mar, ind_i) == geno(mar, ind_j))
                        n_matches++;
                }
            }
            result(ind_i,ind_j) = n_matches;
            result(ind_j,ind_i) = n_typed;
        }
    }

    return result;
}

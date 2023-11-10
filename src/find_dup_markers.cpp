// find subsets of markers with identical genotypes

#include "find_dup_markers.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export(".find_dup_markers_notexact")]]
IntegerVector find_dup_markers_notexact(const IntegerMatrix& Geno, // matrix of genotypes, individuals x markers
                                        const IntegerVector& order, // vector indicating order to be considered, most data to least
                                        const IntegerVector& markerloc, //
                                        const bool adjacent_only)   // if true, consider only adjacent markers
{
    const int n_ind = Geno.rows();
    const int n_mar = Geno.cols();
    if(order.size() != n_mar)
        throw std::invalid_argument("length(order) != ncol(Geno)");
    if(markerloc.size() != n_mar)
        throw std::invalid_argument("length(markerloc) != ncol(Geno)");

    IntegerVector result(n_mar);
    for(int i=0; i<n_mar; i++) result[i] = 0;

    for(int i=0; i<n_mar-1; i++) {
        int oi = order[i]-1;
        for(int j=(i+1); j<n_mar; j++) {
            int oj = order[j]-1;

            if(result[oj] != 0 ||
               (adjacent_only && abs(markerloc[oi] - markerloc[oj]) > 1)) {
                /* skip */
            }
            else {
                int flag = 0;
                for(int k=0; k<n_ind; k++) {
                    if((Geno(k,oi)==0 && Geno(k,oj)!=0) ||
                       (Geno(k,oi)!=0 && Geno(k,oj)!=0 && Geno(k,oi) != Geno(k,oj))) {
                        flag = 1;
                        break;
                    }
                }
                if(!flag) { /* it worked */
                    if(result[oi] != 0) result[oj] = result[oi];
                    else result[oj] = oi+1;
                }
            }
        }
    }

    return(result);
}

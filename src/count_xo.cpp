// count number of crossovers

#include "count_xo.h"
#include "cross.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(".count_xo")]]
IntegerVector count_xo(const IntegerMatrix geno, // genotype matrix markers x individuals
                       const String& crosstype,
                       const bool is_X_chr)
{
    const int n_ind = geno.cols();
    const int n_mar = geno.rows();

    QTLCross* cross = QTLCross::Create(crosstype);

    IntegerVector result(n_ind);
    IntegerVector null_cross_info;

    for(int ind=0; ind<n_ind; ind++) {
        int last_g = 0;
        int n_xo = 0;
        for(int mar=0; mar<n_mar; mar++) {
            int g=geno(mar,ind);
            if(IntegerVector::is_na(g) || g==0) continue; // missing value

            if(last_g==0) { // haven't seen one yet
                last_g = geno(mar,ind);
                continue;
            }

            if(g != last_g)
                n_xo += cross->nrec(last_g, g, is_X_chr, false, null_cross_info);

            last_g = g;
        } // end loop over markers

        result[ind] = n_xo;
    } // end loop over individuals

    delete cross;
    return result;
}

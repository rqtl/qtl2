// locate crossovers

#include "locate_xo.h"
#include "cross.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(".locate_xo")]]
List locate_xo(const IntegerMatrix geno, // genotype matrix markers x individuals
               const NumericVector map,
               const String& crosstype,
               const bool is_X_chr)
{
    const int n_ind = geno.cols();
    const int n_mar = geno.rows();
    if(n_mar != map.size())
        throw std::invalid_argument("Different no. markers in geno and map");

    QTLCross* cross = QTLCross::Create(crosstype);

    std::vector< std::vector<double> > result(n_ind);
    IntegerVector null_cross_info;

    for(int ind=0; ind<n_ind; ind++) {
        int last_g = 0;
        double last_pos = 0.0;
        for(int mar=0; mar<n_mar; mar++) {
            int g=geno(mar,ind);
            if(IntegerVector::is_na(g) || g==0) continue; // missing value

            if(last_g==0) { // haven't seen one yet
                last_g = geno(mar,ind);
                last_pos = map[mar];
                continue;
            }

            if(g != last_g) {
                int n_xo = cross->nrec(last_g, g, is_X_chr, false, null_cross_info);
                for(int xo=0; xo<n_xo; xo++)
                    result[ind].push_back((map[mar] + last_pos)/2.0);
            }

            last_g = g;
            last_pos = map[mar];
        } // end loop over markers
    } // end loop over individuals

    delete cross;
    return wrap(result);
}

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

// [[Rcpp::export(".count_xo_3d")]]
IntegerMatrix count_xo_3d(const IntegerVector geno_array, // 3d array of genotypes, markers x individuals x imputations
                          const String& crosstype,
                          const bool is_X_chr)
{
    if(Rf_isNull(geno_array.attr("dim")))
        throw std::invalid_argument("geno_array should be a 3d array but has no dim attribute");
    const IntegerVector& dim = geno_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("geno_array should be 3d array of genotypes");
    const int n_pos = dim[0];
    const int n_ind = dim[1];
    const int n_imp = dim[2];
    const int mat_size = n_pos*n_ind;

    IntegerMatrix result(n_ind, n_imp);

    for(int imp=0, offset=0; imp<n_imp; imp++, offset += mat_size) {
        IntegerMatrix geno_matrix(n_pos, n_ind);
        std::copy(geno_array.begin()+offset, geno_array.begin()+offset+mat_size, geno_matrix.begin());

        result(_,imp) = count_xo(geno_matrix, crosstype, is_X_chr);
    }

    return result;
}

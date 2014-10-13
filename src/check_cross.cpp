// functions to check QTL cross data/information

#include <Rcpp.h>
#include "cross.h"

using namespace Rcpp;

// check if a cross type is supported
// [[Rcpp::export(".crosstype_supported")]]
bool crosstype_supported(const String& crosstype)
{
    QTLCross* cross = QTLCross::Create(crosstype);
    
    return cross->crosstype_supported();
}

// count inconsistencies in marker data
// [[Rcpp::export(".count_invalid_genotypes")]]
IntegerVector count_invalid_genotypes(const String& crosstype,
                                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                      const bool& is_X_chr,
                                      const LogicalVector& is_female,
                                      const IntegerMatrix& cross_info) // columns are individuals
{
    QTLCross* cross = QTLCross::Create(crosstype);

    int n_ind = genotypes.cols();
    int n_mar = genotypes.rows();

    if(is_female.size() != n_ind)
        throw std::range_error("length(is_female) != ncol(genotypes)");
    if(cross_info.cols() != n_ind)
        throw std::range_error("ncols(cross_info) != ncol(genotypes)");

    IntegerVector result(n_ind);

    for(int ind=0; ind<n_ind; ind++) {
        for(int mar=0; mar<n_mar; mar++) // counting valid genotypes
            result[ind] += cross->check_geno(genotypes(mar,ind), true, is_X_chr, 
                                             is_female[ind], cross_info(_, ind));
        result[ind] = n_mar - result[ind];
    }

    return result;
}                                      
                                      

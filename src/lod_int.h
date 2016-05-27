// calculate lod support intervals
#ifndef LOD_INT_H
#define LOD_INT_H

#include <Rcpp.h>

// this is the "plain" version ignoring the possibility of multiple LOD peaks
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
std::vector<int> lod_int_plain(const Rcpp::NumericVector& lod,
                               const double drop);

// now the version to deal with a chromosome with multiple peaks
//
// well, first a version with a given peak
// input is
//     lod       : vector of lod scores
//     peakindex : index (in 0,1,...,n-1) of peak
//     drop      : amount to drop for support interval
//     peakdrop  : amount to drop between peaks
// this really only makes sense if peakdrop > drop
// output is just like lod_int_plain
//
std::vector<int> lod_int_peak(const Rcpp::NumericVector& lod,
                              const double peakindex,
                              const double peakdrop,
                              const double drop);

#endif // LOD_INT_H

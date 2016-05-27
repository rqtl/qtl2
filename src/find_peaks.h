// find peaks in a lod curve for one chromosome
// above threshold, and dropping by at least peakdrop between peaks
#ifndef FIND_PEAKS_H
#define FIND_PEAKS_H

#include <Rcpp.h>

Rcpp::IntegerVector find_peaks(const Rcpp::NumericVector& lod,
                               const double threshold,
                               const double peakdrop);

#endif // FIND_PEAKS_H

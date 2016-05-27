// find peaks in a lod curve for one chromosome
// above threshold, and dropping by at least peakdrop between peaks

#include "find_peaks.h"
#include <Rcpp.h>
using namespace Rcpp;

// input is a vector of LOD scores ordered by position along a chromosome
// output is a vector of indexes (in 1, 2, 3, ..., lod.size()) of peak locations
// [[Rcpp::export(".find_peaks")]]
IntegerVector find_peaks(const NumericVector& lod, const double threshold, const double peakdrop)
{
    IntegerVector peaks;
    const int n = lod.size();
    int n_peaks = 0;

    double last_peak = 0.0;
    double min_since = 0.0;
    for(int i=0; i<n; i++) {
        if(lod[i] < min_since)
            min_since = lod[i];

        if(lod[i] > threshold) { // possible peak
            if(n_peaks==0) { // first one
                peaks.push_back(i+1);
                last_peak = lod[i];
                min_since = lod[i];
                n_peaks++;
            }
            else { // not the first peak
                if(last_peak > min_since + peakdrop && lod[i] > min_since + peakdrop) { // new peak
                    peaks.push_back(i+1);
                    last_peak = lod[i];
                    min_since = lod[i];
                    n_peaks++;
                }
                else if(lod[i] > last_peak) { // move the last peak
                    peaks[n_peaks-1] = i+1;
                    last_peak = lod[i];
                    min_since = lod[i];
                }
            }
        }
    }

    return peaks;
}

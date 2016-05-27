// find peaks in a lod curve for one chromosome
// above threshold, and dropping by at least peakdrop between peaks

#include "find_peaks.h"
#include <Rcpp.h>
#include "lod_int.h"
using namespace Rcpp;

// this version returns a single index for each peak
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a vector of indexes (in 0, 1, 2, ..., lod.size()-1) of peak locations
//
std::vector<int> find_peaks_plain(const NumericVector& lod,
                                  const double threshold,
                                  const double peakdrop)
{
    std::vector<int> peaks;
    const int n = lod.size();
    int n_peaks = 0;

    double last_peak = 0.0;
    double min_since = 0.0;
    for(int i=0; i<n; i++) {
        if(lod[i] < min_since)
            min_since = lod[i];

        if(lod[i] > threshold) { // possible peak
            if(n_peaks==0) { // first one
                peaks.push_back(i);
                last_peak = lod[i];
                min_since = lod[i];
                n_peaks++;
            }
            else { // not the first peak
                if(last_peak > min_since + peakdrop && lod[i] > min_since + peakdrop) { // new peak
                    peaks.push_back(i);
                    last_peak = lod[i];
                    min_since = lod[i];
                    n_peaks++;
                }
                else if(lod[i] > last_peak) { // move the last peak
                    peaks[n_peaks-1] = i;
                    last_peak = lod[i];
                    min_since = lod[i];
                }
            }
        }
    }

    return peaks;
}



// this version deals with ties in the LOD scores
// ...it returns all indexes which jointly achieve the maximum LOD score
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a list of vectors of indexes (in 0, 1, 2, ..., lod.size()-1) of peak locations
//
// [[Rcpp::export(".find_peaks")]]
std::vector< std::vector<int> > find_peaks(const NumericVector& lod,
                                           const double threshold,
                                           const double peakdrop)
{
    const int n = lod.size();

    // first find the peaks
    std::vector<int> main_peaks = find_peaks_plain(lod, threshold, peakdrop);
    const int n_peaks = main_peaks.size();

    // then go left and right and find adjacent positions with common value
    std::vector< std::vector<int> > result;
    for(int i=0; i<n_peaks; i++) {
        std::vector<int> this_peak;
        this_peak.push_back(main_peaks[i]);

        double this_lod = lod[main_peaks[i]];
        for(int left=main_peaks[i]-1; left>=0; left--) {
            if(lod[left] == this_lod)
                this_peak.push_back(left);
            if(lod[left] < this_lod) break;
        }
        for(int right=main_peaks[i]+1; right<n; right++) {
            if(lod[right] == this_lod)
                this_peak.push_back(right);
            if(lod[right] < this_lod) break;
        }

        result.push_back(this_peak);
    }

    return result;
}


// this version returns both peaks and lod intervals
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a list of vectors of indexes (in 0, 1, 2, ..., lod.size()-1)
//     first two values are the left and right endpoints of the interval
//     remaining values are the indexes with the maximum LOD score
//
// [[Rcpp::export(".find_peaks_and_lodint")]]
std::vector< std::vector<int> > find_peaks_and_lodint(const NumericVector& lod,
                                                      const double threshold,
                                                      const double peakdrop,
                                                      const double drop)
{
    if(drop > peakdrop)
        throw std::invalid_argument("Must have drop <= peakdrop");

    // first find the multiple peaks (can use the simple version)
    std::vector<int> main_peaks = find_peaks_plain(lod, threshold, peakdrop);
    const int n_peaks = main_peaks.size();

    std::vector< std::vector<int> > result;

    // then for each peak, find the lod interval
    for(int i=0; i<n_peaks; i++) {
        std::vector<int> lodint = lod_int_peak(lod,
                                               main_peaks[i],
                                               peakdrop,
                                               drop);

        result.push_back(lodint);
    }

    return(result);
}

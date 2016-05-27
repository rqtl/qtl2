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



// like the plain version, but also returning the locations of the valleys in-between
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a vector of two vectors of indexes (in 0, 1, 2, ..., lod.size()-1)
//    - peak locations
//    - valleys between peaks (including 0 and n-1)
//
std::vector< std::vector<int> > find_peaks_valleys(const NumericVector& lod,
                                                   const double threshold,
                                                   const double peakdrop)
{
    std::vector<int> peaks;
    std::vector<int> valleys;
    const int n = lod.size();
    int n_peaks = 0;
    valleys.push_back(0);

    double last_peak = 0.0;
    double min_since = 0.0;
    int min_since_loc = 0;
    for(int i=0; i<n; i++) {
        if(lod[i] < min_since) {
            min_since = lod[i];
            min_since_loc = i;
        }

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
                    valleys.push_back(min_since_loc);
                    min_since_loc = i;
                    last_peak = lod[i];
                    min_since = lod[i];
                    n_peaks++;
                }
                else if(lod[i] > last_peak) { // move the last peak
                    peaks[n_peaks-1] = i;
                    last_peak = lod[i];
                    min_since = lod[i];
                    min_since_loc = i;
                }
            }
        }
    }
    valleys.push_back(n-1);

    std::vector< std::vector<int> > result;
    result.push_back(peaks);
    result.push_back(valleys);

    return result;
}



// this version deals with ties in the LOD scores
// ...it returns all indexes which jointly achieve the maximum LOD score
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a list of vectors of indexes (in 0, 1, 2, ..., lod.size()-1) of peak locations
//
// The R_ version is a wrapper for R
//
// [[Rcpp::export(".find_peaks")]]
List R_find_peaks(const NumericVector &lod,
                  const double threshold,
                  const double peakdrop)
{
    std::vector< std::vector<int> > result = find_peaks(lod, threshold, peakdrop);

    return wrap(result);
}

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
// The R_ version is a wrapper for R
//
// [[Rcpp::export(".find_peaks_and_lodint")]]
List R_find_peaks_and_lodint(const NumericVector &lod,
                             const double threshold,
                             const double peakdrop,
                             const double drop)
{
    std::vector< std::vector<int> > result = find_peaks_and_lodint(lod, threshold, peakdrop, drop);

    return wrap(result);
}


std::vector< std::vector<int> > find_peaks_and_lodint(const NumericVector& lod,
                                                      const double threshold,
                                                      const double peakdrop,
                                                      const double drop)
{
    if(drop > peakdrop)
        throw std::invalid_argument("Must have drop <= peakdrop");

    // first find the multiple peaks (and valleys between them)
    std::vector< std::vector<int> > peaks_valleys = find_peaks_valleys(lod, threshold, peakdrop);
    std::vector<int> main_peaks = peaks_valleys[0];
    std::vector<int> valleys = peaks_valleys[1];
    const int n_peaks = main_peaks.size();

    std::vector< std::vector<int> > result;

    // then for each peak, find the lod interval
    for(int i=0; i<n_peaks; i++) {
        // find the lod interval, only looking between the adjacent valleys
        std::vector<int> lodint =
            lod_int_contained(lod, main_peaks[i], drop, valleys[i], valleys[i+1]);

        result.push_back(lodint);
    }

    return(result);
}

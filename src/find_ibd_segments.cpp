// identify IBD segments in founder genotypes

#include "find_ibd_segments.h"
#include <Rcpp.h>
using namespace Rcpp;

// find_ibd_segments
// For a pair of individuals on a single chromosome:
//   calculate LOD score for each interval for evidence of IBD vs not
//   find set of non-overlapping intervals with the largest LOD scores
//
// Input:
//   g1  Genotypes of individual 1 (values 1/3, no NAs)
//   g2  Genotypes of individual 2 (values 1/3, no NAs)
//   p   Frequency of genotype 1 at each marker
//   error_prob  Probability of error or mutation at a marker
//   (g1, g2, p all have common length, n = number of markers
//
// Output:
//   matrix with n rows and the following 6 columns
//    - left endpoint (just the values 1...n)
//    - right endpoint (as values in 1...n)
//    - LOD score
//    - number of markers
//    - number of mismatches
//    - 0/1 indicating whether to retain (0 = some overlapping interval has larger LOD score)
//
// [[Rcpp::export(".find_ibd_segments")]]
NumericMatrix find_ibd_segments(const IntegerVector& g1,
                                const IntegerVector& g2,
                                const NumericVector& p,
                                const double error_prob)
{
    const unsigned int n = g1.size();
    if(g2.size() != n)
        throw std::invalid_argument("length(g1) != length(g2)");
    if(p.size() != n)
        throw std::invalid_argument("length(g1) != length(p)");

    NumericMatrix result(n, 6);

    const double log10_error_prob = log10(error_prob);

    // LOD scores at each marker for all pairs of strains
    NumericVector marker_lod(n);
    IntegerVector mismatch(n);
    for(unsigned int i=0; i<n; i++) {
        if(g1[i] == g2[i]) {
            mismatch[i] = 0;
            if(g1[i] == 1) {
                marker_lod[i] = log10((1.0 - error_prob)/p[i] + error_prob);
            } else {
                marker_lod[i] = log10((1.0 - error_prob)/(1.0 - p[i]) + error_prob);
            }
        } else {
            mismatch[i] = 1;
            marker_lod[i] = log10_error_prob;
        }
    }

    // for each marker and each strain pair, find interval with that marker as left endpoint that has maximum LOD score
    for(unsigned int i=0; i<n; i++) {
        double max_lod = marker_lod[i];
        int max_right = i;
        int max_mismatches = mismatch[i];
        double last_lod = marker_lod[i];
        int last_mismatches = mismatch[i];

        for(unsigned int j=i+1; j<n; j++) {
            double this_lod = last_lod + marker_lod[j];
            int this_mismatches = last_mismatches + mismatch[j];
            if(this_lod > max_lod) {
                max_lod = this_lod;
                max_right = j;
                max_mismatches = this_mismatches;
            }
            last_lod = this_lod;
            last_mismatches = this_mismatches;
        }

        result(i,0) = (double)(i+1);
        result(i,1) = (double)(max_right + 1);
        result(i,2) = max_lod;
        result(i,3) = (double)(max_right - i + 1);
        result(i,4) = (double)(max_mismatches);
        result(i,5) = 1.0;
    }


    // reduce to non-overlapping intervals
    for(unsigned int i=0; i<n; i++) {
        if(result(i,5) < 0.5) continue; // skip to next
        for(unsigned int j=i+1; j < result(i,1); j++) {
            if(result(j,2) > result(i,2)) {
                result(i,5) = 0.0;
                break;
            } else {
                result(j,5) = 0.0;
            }
        }
    }

    return(result);
}

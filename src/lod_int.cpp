// calculate lod support intervals

#include "lod_int.h"
#include <Rcpp.h>
#include "find_peaks.h"
using namespace Rcpp;

// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 1, 2, 3, ..., lod.size()) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
// this is the "plain" version ignoring possibility of multiple LOD peaks
// [[Rcpp::export(".lod_int_plain")]]
IntegerVector lod_int_plain(const NumericVector& lod, const double drop)
{
    const int n = lod.size();

    // pass through once to find maximum
    double maxlod = 0.0;
    std::vector<int> maxpos;

    for(int i=0; i<n; i++) {
        if(lod[i] > maxlod) {
            maxpos.clear();
            maxpos.push_back(i+1);
            maxlod = lod[i];
        }
        else if(lod[i] == maxlod) {
            maxpos.push_back(i+1);
        }
    }

    int left, right;
    const double lodmdrop = maxlod - drop;
    for(int i=0, j=n-1; i<n; i++, j--) {
        if(lod[j] > lodmdrop)
            left = j+1;
        if(lod[i] > lodmdrop)
            right = i+1;
    }
    left--;
    right++;
    if(left < 1) left = 1;
    if(right > n) right = n;

    const int n_maxpos = maxpos.size();
    IntegerVector result(n_maxpos+2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}

// now the version to deal with a chromosome with multiple peaks
// well, first a version with a given peak
// input is
//     lod       : vector of lod scores
//     peakindex : index (in 1,...,n) of peak
//     drop      : amount to drop for support interval
//     peakdrop  : amount to drop between peaks
// this really only makes sense if peakdrop > drop
// output is just like lod_int_plain
// [[Rcpp::export(".lod_int_peak")]]
IntegerVector lod_int_peak(const NumericVector& lod,
                           const double peakindex, // index in (1, 2, ..., n) where n = lod.size()
                           const double drop,
                           const double peakdrop)
{
    const int n = lod.size();
    if(peakindex < 1 || peakindex > n)
        throw std::range_error("peakindex out of range");

    if(drop > peakdrop)
        throw std::invalid_argument("Must have drop <= peakdrop");

    double maxlod = lod[peakindex-1];
    std::vector<int> maxpos;
    maxpos.push_back(peakindex);

    const double lodmdrop = maxlod - drop;
    const double lodmpeakdrop = maxlod - peakdrop;
    int left, right;

    // going to the right
    for(int i=peakindex; i<n; i++) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i+1);
        }
        if(lod[i] > lodmdrop) right = i+1;
        if(lod[i] <= lodmpeakdrop) break;
    }
    right++;

    // going to the left
    for(int i=peakindex-2; i>=0; i--) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i+1);
        }
        if(lod[i] > lodmdrop) left = i+1;
        if(lod[i] <= lodmpeakdrop) break;
    }
    left--;

    if(left < 1) left = 1;
    if(right > n) right = n;

    const int n_maxpos = maxpos.size();
    IntegerVector result(n_maxpos+2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}

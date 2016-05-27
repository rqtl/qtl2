// calculate lod support intervals

#include "lod_int.h"
#include <Rcpp.h>
using namespace Rcpp;

// this is the "plain" version ignoring the possibility of multiple LOD peaks
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
// [[Rcpp::export(".lod_int_plain")]]
std::vector<int> lod_int_plain(const NumericVector& lod, const double drop)
{
    const int n = lod.size();

    // pass through once to find maximum
    double maxlod = 0.0;
    std::vector<int> maxpos;

    for(int i=0; i<n; i++) {
        if(lod[i] > maxlod) {
            maxpos.clear();
            maxpos.push_back(i);
            maxlod = lod[i];
        }
        else if(lod[i] == maxlod) {
            maxpos.push_back(i);
        }
    }

    int left, right;
    const double lodmdrop = maxlod - drop;
    for(int i=0, j=n-1; i<n; i++, j--) {
        if(lod[j] > lodmdrop)
            left = j;
        if(lod[i] > lodmdrop)
            right = i;
    }
    left--;
    right++;
    if(left < 0) left = 0;
    if(right > n-1) right = n-1;

    const int n_maxpos = maxpos.size();
    std::vector<int> result(n_maxpos + 2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}

// here we know the peak position and we're looking within a contained subinterval (left, right)
//
// input is a vector of LOD scores ordered by position along a chromosome
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
std::vector<int> lod_int_contained(const NumericVector& lod,
                                   const double peakindex, // index in (0,1,2, ..., n-1) where n = lod.size()
                                   const double drop,
                                   const int start,
                                   const int end)
{
    const int n = lod.size();

    if(peakindex < 0 || peakindex > n-1)
        throw std::range_error("peakindex out of range");
    if(start < 0 || start > n-1)
        throw std::range_error("start out of range");
    if(end < 0 || end > n-1)
        throw std::range_error("end out of range");

    double maxlod = lod[peakindex];
    std::vector<int> maxpos;
    maxpos.push_back(peakindex);

    const double lodmdrop = maxlod - drop;
    int left, right;

    // going to the right
    for(int i=peakindex+1; i<=end; i++) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i);
        }
        if(lod[i] > lodmdrop) right = i;
    }
    right++;

    // going to the left
    for(int i=peakindex-1; i>=start; i--) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i);
        }
        if(lod[i] > lodmdrop) left = i;
    }
    left--;

    if(left < 0) left = 0;
    if(right > n-1) right = n-1;

    // add to vector, at the front
    const int n_maxpos = maxpos.size();
    std::vector<int> result(n_maxpos + 2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}



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
std::vector<int> lod_int_peak(const NumericVector& lod,
                              const double peakindex, // index in (0,1,2, ..., n-1) where n = lod.size()
                              const double peakdrop,
                              const double drop)
{
    const int n = lod.size();
    if(peakindex < 1 || peakindex > n)
        throw std::range_error("peakindex out of range");

    if(drop > peakdrop)
        throw std::invalid_argument("Must have drop <= peakdrop");

    double maxlod = lod[peakindex];
    std::vector<int> maxpos;
    maxpos.push_back(peakindex);

    const double lodmdrop = maxlod - drop;
    const double lodmpeakdrop = maxlod - peakdrop;
    int left, right;

    // going to the right
    for(int i=peakindex+1; i<n; i++) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i);
        }
        if(lod[i] > lodmdrop) right = i;
        if(lod[i] <= lodmpeakdrop) break;
    }
    right++;

    // going to the left
    for(int i=peakindex-1; i>=0; i--) {
        if(lod[i] == maxlod) {
            maxpos.push_back(i);
        }
        if(lod[i] > lodmdrop) left = i;
        if(lod[i] <= lodmpeakdrop) break;
    }
    left--;

    if(left < 0) left = 0;
    if(right > n-1) right = n-1;

    // add to vector, at the front
    const int n_maxpos = maxpos.size();
    std::vector<int> result(n_maxpos + 2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}

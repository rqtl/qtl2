// calculate Bayes credible intervals

#include "bayes_int.h"
#include <Rcpp.h>
#include <algorithm>
#include "util.h"
using namespace Rcpp;


// structures used for sorting area = 10^lod * width
struct area {
    double area;
    int index;
};
// function for sorting the area structures by area
struct by_area {
    bool operator()(const area &a, const area &b) {
        return a.area > b.area; // sorts from largest to smallest
    }
};

// this is the "plain" version ignoring the possibility of multiple LOD peaks
//
// input is a vector of LOD scores ordered by position along a chromosome
// plus the corresponding positions
//
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
// The R_ version is a wrapper for R
//
// [[Rcpp::export(".bayes_int_plain")]]
IntegerVector R_bayes_int_plain(const NumericVector& lod,
                                const NumericVector& pos,
                                const double prob)
{
    std::vector<int> result = bayes_int_plain(lod, pos, prob);

    return wrap(result);
}

std::vector<int> bayes_int_plain(const NumericVector& lod,
                                 const NumericVector& pos,
                                 const double prob)
{
    const int n = lod.size();
    if(n < 2)
        throw std::invalid_argument("Need at least 2 lod scores");
    if(pos.size() != n)
        throw std::invalid_argument("lod and pos should have the same length");

    // calculate interval widths (pos had better be sorted)
    // actually, calculate log(width)
    NumericVector lwidth(n);
    lwidth[0] = log(pos[1] - pos[0]);
    for(int i=1; i<n-1; i++) lwidth[i] = log((pos[i+1] - pos[i-1])/2.0);
    lwidth[n-1] = log(pos[n-1] - pos[n-2]);

    //  10^LOD * interval_width
    std::vector<area> areas(n);
    for(int i=0; i<n; i++) {
        areas[i].area = (lod[i] + lwidth[i])*log(10.0);
        areas[i].index = i;
    }

    // total area
    double total = areas[0].area;
    for(int i=1; i<n; i++)
        total = addlog(total, areas[i].area);

    // sort values from highest to lowest
    std::sort( areas.begin(), areas.end(), by_area());

    int left = n-1;
    int right = 0;
    double cumsum = 0.0;
    for(int i=0; i<n; i++) {
        cumsum += exp(areas[i].area - total);

        int index = areas[i].index;
        if(index < left) left = index;
        if(index > right) right = index;

        if(cumsum >= prob) break; // passed mark
    }

    std::vector<int> result(2);
    result[0] = left;
    result[1] = right;

    return result;
}


// here we know the peak position and we're looking within a contained subinterval (left, right)
//
// input is a vector of LOD scores ordered by position along a chromosome plus
// plus the actual locations
//
// output is a pair of indexes (in 0, 1, 2, ..., lod.size()-1) with endpoints of the interval
// followed by all the indexes where LOD == maximum
//
std::vector<int> bayes_int_contained(const NumericVector& lod,
                                     const NumericVector& pos,
                                     const double peakindex, // index in (0,1,2, ..., n-1) where n = lod.size()
                                     const double prob,
                                     const int start,
                                     const int end)
{
    const int n = lod.size();
    if(n < 2)
        throw std::invalid_argument("Need at least 2 lod scores");
    if(pos.size() != n)
        throw std::invalid_argument("lod and pos should have the same length");

    if(peakindex < 0 || peakindex > n-1)
        throw std::range_error("peakindex out of range");
    if(start < 0 || start > n-1)
        throw std::range_error("start out of range");
    if(end < 0 || end > n-1)
        throw std::range_error("end out of range");
    if(start > end)
        throw std::range_error("should have start <= end");

    double maxlod = lod[peakindex];
    std::vector<int> maxpos;
    maxpos.push_back(peakindex);

    // calculate interval widths (pos had better be sorted)
    // actually, calculate log(width)
    NumericVector lwidth(n);
    lwidth[0] = log(pos[1] - pos[0]);
    for(int i=1; i<n-1; i++) lwidth[i] = log((pos[i+1] - pos[i-1])/2.0);
    lwidth[n-1] = log(pos[n-1] - pos[n-2]);

    const int n_used = end - start + 1;

    //  10^LOD * interval_width
    std::vector<area> areas(n_used);
    for(int i=0; i<n_used; i++) {
        if(start+i != peakindex && lod[start+i] == maxlod)
            maxpos.push_back(peakindex); // finding multiple locations sharing the max LOD
        areas[i].area = (lod[start+i] + lwidth[start+i])*log(10.0);
        areas[i].index = start + i;
    }

    // total area
    double total = areas[0].area;
    for(int i=1; i<n_used; i++)
        total = addlog(total, areas[i].area);

    // sort values from highest to lowest
    std::sort( areas.begin(), areas.end(), by_area());

    int left = peakindex;
    int right = peakindex;
    double cumsum = 0.0;
    for(int i=0; i<n_used; i++) {
        cumsum += exp(areas[i].area - total);

        int index = areas[i].index;
        if(index < left) left = index;
        if(index > right) right = index;

        if(cumsum >= prob) break; // passed mark
    }

    const int n_maxpos = maxpos.size();
    std::vector<int> result(n_maxpos + 2);
    result[0] = left;
    result[1] = right;
    for(int i=0; i<n_maxpos; i++) result[i+2] = maxpos[i];

    return(result);
}

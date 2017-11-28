// Utilities to interpolate between genetic and physical maps

#include "interpolate_maps.h"
#include <exception>
#include <Rcpp.h>
using namespace Rcpp;


// find interval in map that contains pos
// [-1 if to left, map.size()-1 if to right]
// map should be sorted
int find_interval(const double pos, const NumericVector& map)
{
    const int n_map = map.size();

    int result = -1;
    for(int i=0; i<n_map; i++) {
        if(map[i] > pos) return result;
        ++result;
    }
    return result;
}

// for positions relative to oldmap, interpolate to get positions relative to newmap
// [[Rcpp::export]]
NumericVector interpolate_map(const NumericVector& oldpos, const NumericVector& oldmap,
                              const NumericVector& newmap)
{
    const int n_pos = oldpos.size();
    const int n_map = oldmap.size();

    if(newmap.size() != n_map)
        throw std::invalid_argument("length(oldmap) != length(newmap)");

    NumericVector result(n_pos);

    const double oldlength = oldmap[n_map-1] - oldmap[0];
    const double newlength = newmap[n_map-1] - newmap[0];

    for(int i=0; i<n_pos; i++) {
        int interval = find_interval(oldpos[i], oldmap);
        if(interval < 0) {
            if(oldlength == 0.0)
                throw std::invalid_argument("all positions in oldmap coincide");
            result[i] = newmap[0] - (oldmap[0] - oldpos[i])*newlength/oldlength;
        }
        else if(interval == n_map-1) {
            if(oldlength == 0.0)
                throw std::invalid_argument("all positions in oldmap coincide");
            result[i] = newmap[n_map-1] + (oldpos[i] - oldmap[n_map-1])*newlength/oldlength;
        }
        else {
            result[i] = newmap[interval] + (oldpos[i] - oldmap[interval]) *
                (newmap[interval+1] - newmap[interval]) / (oldmap[interval+1] - oldmap[interval]);
        }
    }

    return result;
}

// apply find_interval() to each of a vector of positions
//
// result has two columns and length(pos) rows
//     1st column contains the intervals containing pos
//     2nd column contains 0/1 indicators of whether pos matches left endpoint
//         (to within tolerance tol)
//
// [[Rcpp::export]]
IntegerMatrix find_intervals(const NumericVector& pos,
                             const NumericVector& map,
                             const double tol=1e-8)
{
    const int n_pos = pos.size();
    const int n_map = map.size();
    IntegerMatrix result(n_pos,2);

    for(int i=0; i<n_pos; i++) {
        int interval = find_interval(pos[i], map);
        result(i,0) = interval;

        if(interval < 0 || interval > n_map-1 ||
           fabs(map[interval] - pos[i]) > tol)
            result(i,1) = 0;
        else result(i,1) = 1;
    }

    colnames(result) = CharacterVector::create("interval", "on_map");

    return result;
}

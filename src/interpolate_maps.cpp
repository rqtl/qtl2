// Utilities to interpolate between genetic and physical maps

#include <exception>
#include <Rcpp.h>
using namespace Rcpp;

#include "interpolate_maps.h"

// find interval in map that contains pos
// [-1 if to left, map.size()-1 if to right]
// map should be sorted
int find_interval(const double pos, const NumericVector& map)
{
    const int n_map = map.size();

    int result = -1;
    for(unsigned int i=0; i<n_map; i++) {
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
    int n_pos = oldpos.size();
    int n_map = oldmap.size();

    if(newmap.size() != n_map)
        throw std::invalid_argument("length(oldmap) != length(newmap)");

    NumericVector result(n_pos);

    double oldlength = oldmap[n_map-1] - oldmap[0];
    double newlength = newmap[n_map-1] - newmap[0];

    for(unsigned int i=0; i<n_pos; i++) {
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

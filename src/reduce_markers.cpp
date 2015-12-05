// reduce markers to a more evenly-spaced set

#include "reduce_markers.h"
#include <Rcpp.h>
#include "random.h"

using namespace Rcpp;

// pos = vector of marker positions
// Seek subset of markers such that no two markers are within min_dist
// and sum(weights) is as large as possible
//
// return value is integer vector of marker indices,
// in {1, 2, ..., length(pos)}

// [[Rcpp::export(".reduce_markers")]]
IntegerVector reduce_markers(const NumericVector& pos,
                             const NumericVector& weights,
                             const double min_dist)
{
    const unsigned int n = pos.size();

    if(weights.size() != n)
        throw std::range_error("length(pos) != length(weights)");

    NumericVector total_weights(n);
    IntegerVector prev_marker(n);
    IntegerVector max_to_choose(n);
    IntegerVector path(n);
    unsigned int n_path=0;
    int n_max_to_choose, themax;

    prev_marker[0] = -1;
    total_weights[0] = weights[0];

    for(unsigned int i=1; i<n; i++) {
        if(pos[i] < pos[0] + min_dist) {
            // no markers to left of i that are > min_dist away
            total_weights[i] = weights[i];
            prev_marker[i] = -1;
        }
        else {

            // look for maxima
            n_max_to_choose = 1;
            max_to_choose[0] = 0;
            themax = total_weights[0];
            for(unsigned int j=1; j<i; j++) {

                if(pos[i] < pos[j] + min_dist) break;

                if(total_weights[j] > themax) {
                    n_max_to_choose = 1;
                    max_to_choose[0] = j;
                    themax = total_weights[j];
                }
                else if(total_weights[j] == themax) {
                    max_to_choose[n_max_to_choose] = j;
                    n_max_to_choose++;
                }
            }

            // now choose among the maxima at random
            total_weights[i] = themax + weights[i];
            if(n_max_to_choose == 1) prev_marker[i] = max_to_choose[0];
            else // pick random
                prev_marker[i] = max_to_choose[sample_int(n_max_to_choose)];
        }
    }

    // now find global max
    themax = total_weights[0];
    n_max_to_choose = 1;
    max_to_choose[0] = 0;

    for(unsigned int i=1; i<n; i++) {
        if(total_weights[i] > themax) {
            themax = total_weights[i];
            n_max_to_choose = 1;
            max_to_choose[0] = i;
        }
        else if(total_weights[i] == themax) {
            max_to_choose[n_max_to_choose] = i;
            n_max_to_choose++;
        }
    }

    // right-most marker at global maximum
    if(n_max_to_choose == 1) path[0] = max_to_choose[0];
    else // pick random
        path[0] = max_to_choose[sample_int(n_max_to_choose)];

    n_path=1;

    // trace back
    while(prev_marker[path[n_path-1]] > -1) {
        path[n_path] = prev_marker[path[n_path-1]];
        n_path++;
    }

    // note: result is backwards and has indices starting at 0
    IntegerVector result(n_path);
    for(unsigned int i=0; i<n_path; i++)
        result[i] = path[n_path-1-i] + 1;
    return result;
}

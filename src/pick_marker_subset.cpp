// pick subset of well-spaced markers
//
// More preciesly, find subset of markers that maximizes sum(weights)
// subject to the condition that no two adjacent markers are within
// distance d.

#include "pick_marker_subset.h"
#include "random.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".pick_marker_subset")]]
IntegerVector pick_marker_subset(const NumericVector& pos,      // positions of markers
                                 const double min_d,             // minimum position between markers
                                 const NumericVector& weights)  // weights on the markers
{
    const unsigned int n_pos = pos.size();
    if(n_pos != weights.size())
        throw std::range_error("length(pos) != length(weights)");

    NumericVector total_weights(n_pos);
    IntegerVector prev_marker(n_pos);
    IntegerVector max_to_choose(n_pos);
    IntegerVector path(n_pos);

    unsigned int n_path;
    unsigned int n_max_to_choose;
    double themax;

    /* first location */
    prev_marker[0] = -1;
    total_weights[0] = weights[0];

    for(unsigned int i=1; i<n_pos; i++) {
        if(pos[i] < pos[0] + min_d) {
            /* no markers to left of i that are > min_d away */
            total_weights[i] = weights[i];
            prev_marker[i] = -1;
        }
        else {

            /* look for maxima */
            n_max_to_choose = 1;
            max_to_choose[0] = 0;
            themax = total_weights[0];
            for(unsigned int j=1; j<i; j++) {

                Rcpp::checkUserInterrupt();  // check for ^C from user

                if(pos[i] < pos[j] + min_d) break;

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

            /* now choose among the maxima at random */
            total_weights[i] = themax + weights[i];
            if(n_max_to_choose == 1) prev_marker[i] = max_to_choose[0];
            else /* pick random */
                prev_marker[i] = max_to_choose[random_int(0, n_max_to_choose-1)];
        }
    }

    /* now find global max */
    themax = total_weights[0];
    n_max_to_choose = 1;
    max_to_choose[0] = 0;

    for(unsigned int i=1; i<n_pos; i++) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

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

    /* right-most marker at global maximum */
    if(n_max_to_choose == 1) path[0] = max_to_choose[0];
    else /* pick random */
        path[0] = max_to_choose[random_int(0, n_max_to_choose-1)];

    n_path=1;

    /* trace back */
    while(prev_marker[path[n_path-1]] > -1) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        path[n_path] = prev_marker[path[n_path-1]];
        (n_path)++;
    }

    // the results in "path" are backward and have indexes starting at 0
    IntegerVector result(n_path);
    for(unsigned int i=0; i<n_path; i++)
        result[i] = path[n_path-i-1]+1;

    return result;
}

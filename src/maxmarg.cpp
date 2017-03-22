// identify genotype with maximum marginal probability, arg max Pr(g)

#include "maxmarg.h"
#include <Rcpp.h>
using namespace Rcpp;

// input array as prob[gen,pos,ind]
//     if pr is output of calc_genoprob, need aperm(pr$prob[[1]], c(2,3,1))
// output array as prob[ind, pos]

// [[Rcpp::export(".maxmarg")]]
IntegerMatrix maxmarg(const NumericVector& prob_array, const double minprob)
{
    if(Rf_isNull(prob_array.attr("dim")))
        throw std::invalid_argument("prob_array has no dimension attribute");
    const IntegerVector& dim = prob_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("prob_array should be 3-dimensional array of probabilities");
    const int n_gen = dim[0];
    const int n_pos = dim[1];
    const int n_ind = dim[2];

    IntegerMatrix result(n_ind, n_pos);

    for(int ind=0, offset=0; ind<n_ind; ++ind) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        for(int pos=0; pos<n_pos; pos++, offset += n_gen) {
            double maxp = minprob;
            int state = NA_INTEGER;
            for(int gen=0; gen<n_gen; gen++) {
                double p = prob_array[gen + offset];
                if(p > maxp) {
                    maxp = p;
                    state = gen+1;
                }
            }
            result(ind, pos) = state;
        } // loop over positions
    } // loop over individuals

    return result;
}

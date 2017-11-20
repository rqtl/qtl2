// identify genotype with maximum marginal probability, arg max Pr(g)

#include "maxmarg.h"
#include <Rcpp.h>
#include "random.h"
using namespace Rcpp;

// input array as prob[gen,pos,ind]
//     if pr is output of calc_genoprob, need aperm(pr$prob[[1]], c(2,3,1))
// output array as prob[ind, pos]
//
// states with probabilities within tol are treated as having the same probability

// [[Rcpp::export(".maxmarg")]]
IntegerMatrix maxmarg(const NumericVector& prob_array, const double minprob, const double tol)
{
    if(Rf_isNull(prob_array.attr("dim")))
        throw std::invalid_argument("prob_array should be a 3d array but has no dim attribute");
    const IntegerVector& dim = prob_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("prob_array should be a 3d array of probabilities");
    const int n_gen = dim[0];
    const int n_pos = dim[1];
    const int n_ind = dim[2];

    IntegerMatrix result(n_ind, n_pos);

    for(int ind=0, offset=0; ind<n_ind; ++ind) {
        Rcpp::checkUserInterrupt();  // check for ^C from user

        for(int pos=0; pos<n_pos; pos++, offset += n_gen) {
            double maxp = minprob;
            std::vector<int> states;
            for(int gen=0; gen<n_gen; gen++) {
                double p = prob_array[gen + offset];
                if(fabs(p - maxp) < tol) { // treat the same
                    states.push_back(gen+1);
                }
                else if(p > maxp) {
                    maxp = p;
                    states.clear();
                    states.push_back(gen+1);
                }
            }

            const int n_states = states.size();
            if(n_states==0)
                result(ind, pos) = NA_INTEGER;
            else if(n_states==1)
                result(ind, pos) = states[0];
            else // multiple states with maximum probability; return random choice
                result(ind, pos) = states[sample_int(n_states)];

        } // loop over positions
    } // loop over individuals

    return result;
}

// calculate genetic similarity from genotype probabilities

#include "calc_genetic_sim.h"
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(".calc_genetic_sim")]]
NumericMatrix calc_genetic_sim(const NumericVector& prob_array) // array as n_pos x n_gen x n_ind
{
    if(Rf_isNull(prob_array.attr("dim")))
        throw std::invalid_argument("prob_array has no dimension attribute");
    const IntegerVector& dim = prob_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("prob_array should be 3-dimensional array of probabilities");
    const unsigned int n_ind = dim[2];
    const unsigned int n_pos = dim[0];
    const unsigned int n_gen = dim[1];
    const unsigned int pos_by_gen = n_pos*n_gen;

    NumericMatrix result(n_ind, n_ind);

    for(unsigned int ind_i=0, offset_i=0; ind_i<n_ind; ++ind_i, offset_i += pos_by_gen) {
        result(ind_i,ind_i) = (double)n_pos;

        for(unsigned int ind_j=ind_i+1, offset_j=(ind_i+1)*pos_by_gen; ind_j<n_ind; ind_j++, offset_j += pos_by_gen) {

            double total = 0.0;
            for(unsigned int pos=0; pos<n_pos; pos++) {
                for(unsigned int gen=0; gen<n_gen; gen++) {
                    total += prob_array[pos*n_gen + gen + offset_i] *
                        prob_array[pos*n_gen + gen + offset_j];
                }
            }
            result(ind_i,ind_j) = result(ind_j,ind_i) = total;
        }
    }

    return result;
}

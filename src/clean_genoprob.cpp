// clean genoprobs, setting small values to 0

#include "clean_genoprob.h"
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;

// clean genoprobs, setting small values to 0
// [[Rcpp::export(".clean_genoprob")]]
NumericVector clean_genoprob(const NumericVector& prob_array, // array as n_ind x n_gen x n_pos
                             double value_threshold=1e-6,
                             double column_threshold=0.01)
{
    if(Rf_isNull(prob_array.attr("dim")))
        throw std::invalid_argument("prob_array should be a 3d array but has no dimension attribute");
    const IntegerVector& dim = prob_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("prob_array should be a 3d array of probabilities");
    const int n_ind = dim[0];
    const int n_gen = dim[1];
    const int n_pos = dim[2];

    NumericVector result = clone(prob_array);

    // ensure that we don't set all values in a row to 0
    if(column_threshold > 1.0/(double)n_gen)
        column_threshold = 0.5/(double)n_gen;
    if(value_threshold > 1.0/(double)n_gen)
        value_threshold = 0.5/(double)n_gen;

    for(int pos=0, offset=0; pos<n_pos; pos++) {

        // first look at each genotype column and find max; if < column_threshold, set all values to 0
        for(int gen=0; gen<n_gen; gen++) {
            bool zero_column = true;
            for(int ind=0; ind<n_ind; ind++) {
                if(prob_array[ind + gen*n_ind + pos*n_gen*n_ind] >= column_threshold) {
                    zero_column = false;
                    break;
                }
            }
            if(zero_column) { // biggest value was < column_threshold so zero the column
                for(int ind=0; ind<n_ind; ind++) {
                    result[ind + gen*n_ind + pos*n_gen*n_ind] = 0.0;
                }
            }
        }

        // now look at the individual values

        for(int ind=0; ind<n_ind; ind++) {
            double sum=0.0;
            for(int gen=0; gen<n_gen; gen++) {
                // small values set to 0
                int index = ind + gen*n_ind + pos*n_gen*n_ind;
                if(result[index] < value_threshold) result[offset+gen] = 0.0;

                // get sum so we can rescale to sum to 1
                sum += result[offset+gen];
            }

            for(int gen=0; gen<n_gen; gen++) {
                int index = ind + gen*n_ind + pos*n_gen*n_ind;
                result[offset+gen] /= sum;
            }

        }
    }

    result.attr("dim") = Dimension(n_ind, n_gen, n_pos);

    return result;
}

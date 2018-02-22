// threshold genoprobs, setting small values to 0

#include "threshold_genoprob.h"
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;

// threshold genoprobs, setting small values to 0
// [[Rcpp::export(".threshold_genoprob")]]
NumericVector threshold_genoprob(const NumericVector& prob_array, // array as n_gen x n_ind x n_pos
                                 const double threshold=1e-6)
{
    if(Rf_isNull(prob_array.attr("dim")))
        throw std::invalid_argument("prob_array should be a 3d array but has no dimension attribute");
    const IntegerVector& dim = prob_array.attr("dim");
    if(dim.size() != 3)
        throw std::invalid_argument("prob_array should be a 3d array of probabilities");
    const int n_gen = dim[0];
    const int n_ind = dim[1];
    const int n_pos = dim[2];
    const int ind_by_pos = n_ind*n_pos;

    NumericVector result(n_gen*n_ind*n_pos);

    for(int pos=0, offset=0; pos<n_pos; pos++) {
        for(int ind=0; ind<n_ind; ind++, offset += n_gen) {
            double sum=0.0;
            for(int gen=0; gen<n_gen; gen++) {
                // small values set to 0
                if(prob_array[offset+gen] < threshold) result[offset+gen] = 0.0;
                else result[offset+gen] = prob_array[offset+gen];

                // get sum so we can rescale to sum to 1
                sum += result[offset+gen];
            }

            if(sum < threshold*n_gen) { // everything was small ... set each to 1/n_gen
                for(int gen=0; gen<n_gen; gen++) {
                    result[offset+gen] = 1.0/(double)(n_gen);
                }
            }
            else {
                for(int gen=0; gen<n_gen; gen++) {
                    result[offset+gen] /= sum;
                }
            }

        }
    }

    result.attr("dim") = Dimension(n_gen, n_ind, n_pos);

    return result;
}

// interpolate genotype probabilities

#include "interp_genoprob.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".interp_genoprob_onechr")]]
NumericVector interp_genoprob_onechr(const NumericVector& genoprob,
                                     const NumericVector& map,
                                     const IntegerVector& pos_index)
{
    // get dimensions
    if(Rf_isNull(genoprob.attr("dim")))
        throw std::invalid_argument("genoprob should be a 3d array but has no dim attribute");
    const IntegerVector& d = genoprob.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprob should be a 3d array");
    const int n_ind = d[0];
    const int n_gen = d[1];
    const int matsize = n_ind * n_gen;
    const int n_pos = map.size();
    if(pos_index.size() != n_pos) {
        throw std::invalid_argument("Need length(map) == length(pos_index)");
    }

    NumericVector result(n_ind*n_gen*n_pos);
    result.attr("dim") = Dimension(n_ind, n_gen, n_pos);

    // find position to the left that has genoprobs
    IntegerVector left_index(n_pos);
    int last = -1;
    for(int pos=0; pos<n_pos; pos++) {
        if(pos_index[pos] >= 0) last = pos;
        left_index[pos] = last;
    }

    // find position to the right that has genoprobs
    IntegerVector right_index(n_pos);
    last = -1;
    for(int pos=n_pos-1; pos>=0; pos--) {
        if(pos_index[pos] >= 0) last = pos;
        right_index[pos] = last;
    }

    // copy or interpolate
    for(int pos=0; pos<n_pos; pos++) {
        if(pos_index[pos] >= 0) { // in the old genoprobs
            std::copy(genoprob.begin()+(pos_index[pos]*matsize),
                      genoprob.begin()+((pos_index[pos]+1)*matsize),
                      result.begin()+(pos*matsize));
        }
        else {
            double p,q;
            if(left_index[pos] < 0) { // off end to left
                p = 0.0;
                q = 1.0;
            }
            else if(right_index[pos] < 0) { // off end to right
                p = 1.0;
                q = 0.0;
            }
            else {
                double left_pos =  map[left_index[pos]];
                double right_pos = map[right_index[pos]];
                p = (right_pos - map[pos])/(right_pos - left_pos);
                q = (map[pos] - left_pos)/(right_pos - left_pos);
            }

            for(int ind=0; ind<n_ind; ind++) {
                for(int gen=0; gen<n_gen; gen++) {
                    result[ind + gen*n_ind + pos*matsize] = 0.0;
                    if(p > 0)
                        result[ind + gen*n_ind + pos*matsize] +=
                            (p*genoprob[ind + gen*n_ind + pos_index[left_index[pos]]*matsize]);
                    if(q > 0)
                        result[ind + gen*n_ind + pos*matsize] +=
                            (q*genoprob[ind + gen*n_ind + pos_index[right_index[pos]]*matsize]);
                }
            }
        }
    }

    return result;
}

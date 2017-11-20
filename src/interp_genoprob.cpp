// interpolate genotype probabilities

#include "interp_genoprob.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".interp_genoprob_onechr")]]
NumericVector interp_genoprob_onechr(const NumericVector& genoprob,
                                     const NumericVector& map,
                                     const LogicalVector& is_new_pos)
{
    // get dimensions
    const IntegerVector& d = genoprob.attr("dim");
    if(d.size() != 3)
        throw std::invalid_argument("genoprob shoudl be a 3d array");
    const int n_ind = d[0];
    const int n_gen = d[1];
    const int n_mar = d[2];
    const int matsize = n_ind * n_gen;
    const int n_pos = map.size();
    if(is_new_pos.size() != n_pos) {
        throw std::invalid_argument("Need length(map) == length(is_new_pos)");
    }

    NumericVector result(n_ind*n_gen*n_pos);
    result.attr("dim") = Dimension(n_ind, n_gen, n_pos);

    // find position to the right that has genoprobs
    NumericVector right_index(n_pos), right_map(n_pos);
    int last_oldindex=n_mar;
    int last_newindex=n_pos;
    for(int pos=n_pos-1; pos>=0; pos--) {
        if(!is_new_pos[pos]) {
            last_oldindex--;
            last_newindex = pos;
        }
        right_index[pos] = last_oldindex;
        right_map[pos] = map[last_newindex];
    }

    // copy or interpolate
    int left_oldindex = -1;
    int left_newindex = -1;
    for(int pos=0; pos<n_pos; pos++) {
        if(!is_new_pos[pos]) { // in the old genoprobs
            left_oldindex++;
            left_newindex = pos;
            std::copy(genoprob.begin()+(left_oldindex*matsize),
                      genoprob.begin()+((left_oldindex+1)*matsize),
                      result.begin()+(pos*matsize));
        }
        else {
            double leftpos = map[left_oldindex];
            double rightpos = right_map[pos];
            double p, q;
            if(leftpos == rightpos) { p = q = 0.5; }
            else if(right_index[pos] == n_mar) { // off the end to the right
                q = 1.0; // weight on the left value
                p = 0.0; // weight on the right value
            }
            else if(left_oldindex == -1) { // off the end to the left
                p = 1.0; // weight on the right value
                q = 0.0; // weight on the left value
            }
            else {
                p = (map[pos] - leftpos) / (rightpos - leftpos); // weight on the right value
                q = 1.0 - p; // weight on the left value
            }

            for(int ind=0; ind<n_ind; ind++) {
                for(int gen=0; gen<n_ind; gen++) {
                    result[ind + gen*n_ind + pos*matsize] = 0.0;
                    if(q > 0)
                        result[ind + gen*n_ind + pos*matsize] +=
                            (q*genoprob[ind + gen*n_ind + left_oldindex*matsize]);
                    if(p > 0)
                        result[ind + gen*n_ind + pos*matsize] +=
                            (p*genoprob[ind + gen*n_ind + right_index[pos]*matsize]);
                }
            }
        }
    }

    return result;
}

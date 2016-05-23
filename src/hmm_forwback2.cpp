// forward-backward equations for HMM
// (this version assumes constant is_female and cross_info and pre-calcs the step and emit matrices)

#include "hmm_forwback2.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm_util.h"

// forward equations
NumericMatrix forwardEquations2(const IntegerVector& genotypes,
                                const Rcpp::NumericVector& init_vector,
                                const std::vector<Rcpp::NumericMatrix>& emit_matrix,
                                const std::vector<Rcpp::NumericMatrix>& step_matrix,
                                const IntegerVector& marker_index,
                                const IntegerVector& poss_gen)
{
    int n_pos = marker_index.size();

    // possible genotypes for this chromosome and individual
    int n_gen = poss_gen.size();

    // to contain ln Pr(G_i = g | marker data)
    NumericMatrix alpha(n_gen, n_pos);

    // initialize alphas
    for(int i=0; i<n_gen; i++) {
        int g = poss_gen[i];
        alpha(i,0) = init_vector[i];
        if(marker_index[0] >= 0)
            alpha(i,0) += emit_matrix[marker_index[0]](genotypes[marker_index[0]], i);
    }

    for(int pos=1; pos<n_pos; pos++) {
        for(int ir=0; ir<n_gen; ir++) {
            alpha(ir,pos) = alpha(0, pos-1) + step_matrix[pos-1](0, ir);

            for(int il=1; il<n_gen; il++)
                alpha(ir,pos) = addlog(alpha(ir,pos), alpha(il,pos-1) + step_matrix[pos-1](il, ir));

            if(marker_index[pos]>=0)
                alpha(ir,pos) += emit_matrix[marker_index[pos]](genotypes[marker_index[pos]], ir);
        }
    }

    return alpha;
}



// backward Equations
NumericMatrix backwardEquations2(const IntegerVector& genotypes,
                                 const Rcpp::NumericVector& init_vector,
                                 const std::vector<Rcpp::NumericMatrix>& emit_matrix,
                                 const std::vector<Rcpp::NumericMatrix>& step_matrix,
                                 const IntegerVector& marker_index,
                                 const IntegerVector& poss_gen)
{
    int n_pos = marker_index.size();

    // possible genotypes for this chromosome and individual
    int n_gen = poss_gen.size();

    // to contain ln Pr(G_i = g | marker data)
    NumericMatrix beta(n_gen, n_pos);

    // backward equations
    for(int pos = n_pos-2; pos >= 0; pos--) {
        for(int il=0; il<n_gen; il++) {
            for(int ir=0; ir<n_gen; ir++) {
                double to_add = beta(ir,pos+1) + step_matrix[pos](il, ir);
                if(marker_index[pos+1] >=0)
                    to_add += emit_matrix[marker_index[pos+1]](genotypes[marker_index[pos+1]], ir);

                if(ir==0) beta(il,pos) = to_add;
                else beta(il,pos) = addlog(beta(il,pos), to_add);
            }
        }
    }

    return beta;
}

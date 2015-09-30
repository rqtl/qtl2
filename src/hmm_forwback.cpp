// forward-backward equations for HMM

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm_util.h"
#include "hmm_forwback.h"

// forward equations
NumericMatrix forwardEquations(QTLCross* cross,
                               const IntegerVector& genotypes,
                               const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                               const bool is_X_chr,
                               const bool is_female,
                               const IntegerVector& cross_info,
                               const NumericVector& rec_frac,
                               const IntegerVector& marker_index,
                               const double error_prob,
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
        alpha(i,0) = cross->init(g, is_X_chr, is_female, cross_info);
        if(marker_index[0] >= 0)
            alpha(i,0) += cross->emit(genotypes[marker_index[0]], g, error_prob,
                                      founder_geno(_, marker_index[0]), is_X_chr, is_female, cross_info);
    }

    for(int pos=1; pos<n_pos; pos++) {
        for(int ir=0; ir<n_gen; ir++) {
            alpha(ir,pos) = alpha(0, pos-1) + cross->step(poss_gen[0], poss_gen[ir], rec_frac[pos-1],
                                                          is_X_chr, is_female, cross_info);

            for(int il=1; il<n_gen; il++)
                alpha(ir,pos) = addlog(alpha(ir,pos), alpha(il,pos-1) +
                                       cross->step(poss_gen[il], poss_gen[ir], rec_frac[pos-1],
                                                   is_X_chr, is_female, cross_info));

            if(marker_index[pos]>=0)
                alpha(ir,pos) += cross->emit(genotypes[marker_index[pos]], poss_gen[ir], error_prob,
                                             founder_geno(_, marker_index[pos]), is_X_chr, is_female, cross_info);
        }
    }

    return alpha;
}



// backward Equations
NumericMatrix backwardEquations(QTLCross* cross,
                                const IntegerVector& genotypes,
                                const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                const bool is_X_chr,
                                const bool is_female,
                                const IntegerVector& cross_info,
                                const NumericVector& rec_frac,
                                const IntegerVector& marker_index,
                                const double error_prob,
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
                double to_add = beta(ir,pos+1) + cross->step(poss_gen[il], poss_gen[ir], rec_frac[pos],
                                                              is_X_chr, is_female, cross_info);
                if(marker_index[pos+1] >=0)
                    to_add += cross->emit(genotypes[marker_index[pos+1]], poss_gen[ir], error_prob,
                                          founder_geno(_, marker_index[pos+1]), is_X_chr, is_female, cross_info);

                if(ir==0) beta(il,pos) = to_add;
                else beta(il,pos) = addlog(beta(il,pos), to_add);
            }
        }
    }

    return beta;
}

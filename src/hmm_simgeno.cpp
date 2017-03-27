// simulate genotypes given observed marker data

#include "hmm_simgeno.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm_util.h"
#include "hmm_forwback.h"
#include "random.h"

// simulate genotypes given observed marker data
// [[Rcpp::export(".sim_geno")]]
IntegerVector sim_geno(const String& crosstype,
                       const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                       const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                       const bool is_X_chr,
                       const LogicalVector& is_female, // length n_ind
                       const IntegerMatrix& cross_info, // columns are individuals
                       const NumericVector& rec_frac,   // length nrow(genotypes)-1
                       const IntegerVector& marker_index, // length nrow(genotypes)
                       const double error_prob,
                       const int n_draws) // number of imputations
{
    const int n_ind = genotypes.cols();
    const int n_pos = marker_index.size();
    const int n_mar = genotypes.rows();

    QTLCross* cross = QTLCross::Create(crosstype);

    // check inputs
    if(is_female.size() != n_ind)
        throw std::range_error("length(is_female) != ncol(genotypes)");
    if(cross_info.cols() != n_ind)
        throw std::range_error("ncols(cross_info) != ncol(genotypes)");
    if(rec_frac.size() != n_pos-1)
        throw std::range_error("length(rec_frac) != length(marker_index)-1");

    if(error_prob < 0.0 || error_prob > 1.0)
        throw std::range_error("error_prob out of range");

    for(int i=0; i<rec_frac.size(); i++) {
        if(rec_frac[i] < 0 || rec_frac[i] > 0.5)
            throw std::range_error("rec_frac must be >= 0 and <= 0.5");
    }
    if(!cross->check_founder_geno_size(founder_geno, n_mar))
        throw std::range_error("founder_geno is not the right size");
    // end of checks

    const int mat_size = n_pos*n_draws;
    IntegerVector draws(mat_size*n_ind); // output object

    for(int ind=0; ind<n_ind; ind++) {

        Rcpp::checkUserInterrupt();  // check for ^C from user

        // possible genotypes for this individual
        IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
        const int n_poss_gen = poss_gen.size();
        NumericVector probs(n_poss_gen);

        // backward equations
        NumericMatrix beta = backwardEquations(cross, genotypes(_,ind), founder_geno, is_X_chr, is_female[ind],
                                               cross_info(_,ind), rec_frac, marker_index, error_prob,
                                               poss_gen);

        // simulate genotypes
        for(int draw=0; draw<n_draws; draw++) {
            // first draw
            // calculate first prob (on log scale)
            probs[0] = cross->init(poss_gen[0], is_X_chr, is_female[ind], cross_info(_,ind)) + beta(0,0);
            if(marker_index[0] >= 0)
                probs[0] += cross->emit(genotypes(marker_index[0],ind), poss_gen[0], error_prob,
                                        founder_geno(_, marker_index[0]), is_X_chr, is_female[ind], cross_info(_,ind));
            double sumprobs = probs[0]; // to contain log(sum(probs))

            // calculate rest of probs
            for(int g=1; g<n_poss_gen; g++) {
                probs[g] = cross->init(poss_gen[g], is_X_chr, is_female[ind], cross_info(_,ind)) + beta(g,0);
                if(marker_index[0] >= 0)
                    probs[g] += cross->emit(genotypes(marker_index[0],ind), poss_gen[g], error_prob,
                                            founder_geno(_, marker_index[0]), is_X_chr, is_female[ind], cross_info(_,ind));
                sumprobs = addlog(sumprobs, probs[g]);
            }

            // re-scale probs
            for(int g=0; g<n_poss_gen; g++)
                probs[g] = exp(probs[g] - sumprobs);

            // make draw, returns a value from 1, 2, ..., n_poss_gen
            int curgeno = sample_int(probs);
            draws[draw*n_pos + ind*mat_size] = poss_gen[curgeno];

            // move along chromosome
            for(int pos=1; pos<n_pos; pos++) {

                // calculate probs
                for(int g=0; g<n_poss_gen; g++) {
                    probs[g] = cross->step(poss_gen[curgeno], poss_gen[g], rec_frac[pos-1],
                                           is_X_chr, is_female[ind], cross_info(_,ind)) +
                        beta(g,pos) - beta(curgeno, pos-1);
                    if(marker_index[pos] >= 0)
                        probs[g] += cross->emit(genotypes(marker_index[pos],ind), poss_gen[g], error_prob,
                                                founder_geno(_, marker_index[pos]), is_X_chr, is_female[ind], cross_info(_,ind));
                    probs[g] = exp(probs[g]);
                }

                // make draw
                curgeno = sample_int(probs);

                draws[pos + draw*n_pos + ind*mat_size] = poss_gen[curgeno];

            } // loop over positions
        } // loop over draws
    } // loop over individuals

    draws.attr("dim") = Dimension(n_pos, n_draws, n_ind);
    delete cross;
    return draws;
}

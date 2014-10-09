// main HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm.h"

// re-estimate inter-marker recombination fractions
// [[Rcpp::export(".est_map")]]
NumericVector est_map(String crosstype,
                      IntegerMatrix genotypes,
                      bool is_X_chr,
                      LogicalVector is_female,
                      IntegerMatrix cross_info,
                      NumericVector rec_frac,
                      double error_prob,
                      int max_iterations,
                      double tol,
                      bool verbose)
{
    int n_ind = genotypes.cols();
    int n_mar = genotypes.rows();
    int n_rf = n_mar-1;

    // check inputs
    if(is_female.size() != n_ind)
        throw std::range_error("length(is_female) != ncol(genotypes)");
    if(cross_info.cols() != n_ind)
        throw std::range_error("ncols(cross_info) != ncol(genotypes)");
    if(rec_frac.size() != n_rf)
        throw std::range_error("length(rec_frac) != nrow(genotypes)-1");

    if(error_prob < 0.0 || error_prob > 1.0)
        throw std::range_error("error_prob out of range");

    for(int i=0; i<rec_frac.size(); i++) {
        if(rec_frac[i] < 0 || rec_frac[i] > 0.5)
            throw std::range_error("rec_frac must be >= 0 and <= 0.5");
    }

    if(max_iterations < 0)
      throw std::range_error("max_iterations should be >= 0");
    if(tol < 0)
      throw std::range_error("tol >= 0");
    // end of checks

    QTLCross* cross = QTLCross::Create(crosstype);
    if(cross->type != cross->phase_known_type) // get phase-known version of cross
        cross = QTLCross::Create(cross->phase_known_type);
    
    NumericVector cur_rec_frac(n_rf);
    NumericVector prev_rec_frac(clone(rec_frac));

    // marker index for forward/backward equations
    IntegerVector marker_index(n_mar);
    for(int i=0; i<n_mar; i++) marker_index[i] = i;

    // 3-d array to contain sum(gamma(il,ir)) for each interval
    int n_gen = cross->ngen(is_X_chr);
    int n_gen_sq = n_gen*n_gen;
    NumericVector full_gamma(n_gen_sq*n_rf);

    for(int it=0; it<max_iterations; it++) {

        // zero the full_gamma array
        for(int i=0; i<n_gen_sq*n_rf; i++) full_gamma[i] = 0.0;

        for(int ind=0; ind < n_ind; ind++) {

            // possible genotypes for this individual
            IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
            int n_poss_gen = poss_gen.size();

            // forward and backward equations
            NumericMatrix alpha = forwardEquations(cross, genotypes(_,ind), is_X_chr, is_female[ind],
                                                   cross_info(_,ind), prev_rec_frac, marker_index, error_prob,
                                                   poss_gen);
            NumericMatrix beta = backwardEquations(cross, genotypes(_,ind), is_X_chr, is_female[ind],
                                                   cross_info(_,ind), prev_rec_frac, marker_index, error_prob,
                                                   poss_gen);

            for(int pos=0; pos<n_rf; pos++) {
                // calculate gamma = log Pr(v1, v2, O)
                NumericMatrix gamma(alpha.rows(), alpha.cols());
                double sum_gamma;
                bool sum_gamma_undef = true;
                for(int il=0; il<n_poss_gen; il++) {
                    for(int ir=0; ir<n_poss_gen; ir++) {
                        gamma(il,ir) = alpha(il,pos) + beta(ir,pos+1) +
                            cross->emit(genotypes(pos+1,ind), poss_gen[ir], error_prob, is_X_chr, is_female[ind], cross_info(_,ind)) +
                            cross->step(poss_gen[il], poss_gen[ir], prev_rec_frac[pos], is_X_chr, is_female[ind], cross_info(_,ind));

                        if(sum_gamma_undef) {
                            sum_gamma_undef = false;
                            sum_gamma = gamma(il,ir);
                        }
                        else {
                            sum_gamma = addlog(sum_gamma, gamma(il,ir));
                        }
                    }
                }
                
                // add to full_gamma
                for(int il=0; il<n_poss_gen; il++) {
                    int gl = poss_gen[il]-1;
                    for(int ir=0; ir<n_poss_gen; ir++) {
                        int gr = poss_gen[ir]-1;
                        full_gamma[n_gen_sq*pos + gr*n_gen + gl] += exp(gamma(il,ir) - sum_gamma);
                    }
                }
            } // loop over marker intervals
            
        } // loop over individuals

        // re-estimate rec'n fractions
        for(int pos=0; pos < n_rf; pos++) {
            NumericMatrix sub_gamma(n_gen, n_gen);
            std::copy(full_gamma.begin()+n_gen_sq*pos, full_gamma.begin()+n_gen_sq*(pos+1), sub_gamma.begin());
            cur_rec_frac[pos] = cross->est_rec_frac(sub_gamma, is_X_chr);
        }

        if(verbose) {
            double maxdif = max(abs(prev_rec_frac - cur_rec_frac));
            Rprintf("%4d %.12f\n", it, maxdif);
        }

        // check convergence
        bool converged = true;
        for(int pos=0; pos<n_rf; pos++) {
            if(fabs(prev_rec_frac[pos] - cur_rec_frac[pos]) > tol*(cur_rec_frac[pos]+tol*100.0)) {
                converged = false;
                break;
            }
        }

        if(converged) break;

        prev_rec_frac = clone(cur_rec_frac);
    } // end loop over iterations

    // calculate log likelihood
    double loglik = 0.0;
    for(int ind=0; ind<n_ind; ind++) {
        double curloglik;

        // possible genotypes for this individual
        IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
        int n_poss_gen = poss_gen.size();

        // forward and backward equations
        NumericMatrix alpha = forwardEquations(cross, genotypes(_,ind), is_X_chr, is_female[ind],
                                               cross_info(_,ind), cur_rec_frac, marker_index, error_prob,
                                               poss_gen);
        
        bool curloglik_undef = true;
        for(int i=0; i<n_poss_gen; i++) {
            if(curloglik_undef) {
                curloglik_undef = false;
                curloglik = alpha(i,n_rf);
            }
            else {
                curloglik = addlog(curloglik, alpha(i, n_rf));
            }
        }
        loglik += curloglik;
    }

    if(verbose) {
        Rprintf("loglik = %.3f\n", loglik);
    }

    cur_rec_frac.attr("loglik") = loglik;
    return cur_rec_frac;
}

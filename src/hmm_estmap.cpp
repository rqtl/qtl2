// re-estimate inter-marker recombination fractions

#include "hmm_estmap.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm_util.h"
#include "hmm_forwback.h"

// re-estimate inter-marker recombination fractions
// [[Rcpp::export(".est_map")]]
NumericVector est_map(const String& crosstype,
                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                      const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                      const bool is_X_chr,
                      const LogicalVector& is_female,
                      const IntegerMatrix& cross_info,
                      const NumericVector& rec_frac,
                      const double error_prob,
                      const int max_iterations,
                      const double tol,
                      const bool verbose)
{
    int n_ind = genotypes.cols();
    int n_mar = genotypes.rows();
    int n_rf = n_mar-1;

    QTLCross* cross_pu = QTLCross::Create(crosstype);
    QTLCross* cross;
    if(cross_pu->crosstype != cross_pu->phase_known_crosstype) // get phase-known version of cross
        cross = QTLCross::Create(cross_pu->phase_known_crosstype);
    else cross = cross_pu;

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

    if(!cross->check_founder_geno_size(founder_geno, n_mar))
        throw std::range_error("founder_geno is not the right size");
    // end of checks

    NumericVector cur_rec_frac(n_rf);
    NumericVector prev_rec_frac(clone(rec_frac));

    // marker index for forward/backward equations
    IntegerVector marker_index(n_mar);
    for(int i=0; i<n_mar; i++) marker_index[i] = i;

    // 3-d array to contain sum(gamma(il,ir)) for each interval
    int n_gen = cross->ngen(is_X_chr);
    int n_gen_sq = n_gen*n_gen;
    int n_gen_sq_times_n_ind = n_gen_sq * n_ind;
    NumericVector full_gamma(n_gen_sq_times_n_ind * n_rf);

    for(int it=0; it<max_iterations; it++) {

        // zero the full_gamma array
        full_gamma.fill(0.0);

        for(int ind=0; ind < n_ind; ind++) {

            Rcpp::checkUserInterrupt();  // check for ^C from user

            // possible genotypes for this individual
            IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
            int n_poss_gen = poss_gen.size();

            // forward and backward equations
            NumericMatrix alpha = forwardEquations(cross, genotypes(_,ind), founder_geno, is_X_chr, is_female[ind],
                                                   cross_info(_,ind), prev_rec_frac, marker_index, error_prob,
                                                   poss_gen);
            NumericMatrix beta = backwardEquations(cross, genotypes(_,ind), founder_geno, is_X_chr, is_female[ind],
                                                   cross_info(_,ind), prev_rec_frac, marker_index, error_prob,
                                                   poss_gen);

            for(int pos=0; pos<n_rf; pos++) {
                // calculate gamma = log Pr(v1, v2, O)
                NumericMatrix gamma(n_poss_gen, n_poss_gen);
                double sum_gamma=0.0;
                bool sum_gamma_undef = true;
                for(int ir=0; ir<n_poss_gen; ir++) {
                    for(int il=0; il<n_poss_gen; il++) {
                        gamma(il,ir) = alpha(il,pos) + beta(ir,pos+1) +
                            cross->emit(genotypes(pos+1,ind), poss_gen[ir], error_prob,
                                        founder_geno(_,pos+1), is_X_chr, is_female[ind], cross_info(_,ind)) +
                            cross->step(poss_gen[il], poss_gen[ir], prev_rec_frac[pos],
                                        is_X_chr, is_female[ind], cross_info(_,ind));

                        if(sum_gamma_undef) {
                            sum_gamma_undef = false;
                            sum_gamma = gamma(il,ir);
                        }
                        else {
                            sum_gamma = addlog(sum_gamma, gamma(il,ir));
                        }
                    }
                }

                // add to full_gamma array of dim n_rf x n_ind x n_gen x n_gen
                const int offset = n_gen_sq_times_n_ind*pos + n_gen_sq*ind;
                for(int ir=0; ir<n_poss_gen; ir++) {
                    int gr_by_n_gen = (poss_gen[ir]-1)*n_gen;
                    for(int il=0; il<n_poss_gen; il++) {
                        int gl = poss_gen[il]-1;
                        full_gamma[offset + gr_by_n_gen + gl] += exp(gamma(il,ir) - sum_gamma);
                    }
                }
            } // loop over marker intervals

        } // loop over individuals

        // re-estimate rec'n fractions
        for(int pos=0; pos < n_rf; pos++) {
            // pull out the part for that position
            NumericVector sub_gamma(n_gen_sq_times_n_ind);
            std::copy(full_gamma.begin()+(n_gen_sq_times_n_ind*pos),
                      full_gamma.begin()+(n_gen_sq_times_n_ind*(pos+1)),
                      sub_gamma.begin());
            cur_rec_frac[pos] = cross->est_rec_frac(sub_gamma, is_X_chr, cross_info, n_gen);
        }

        if(verbose) {
            double maxdif = max(abs(prev_rec_frac - cur_rec_frac));
            Rprintf("%4d %.12f\n", it+1, maxdif);
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
        double curloglik=0.0;

        Rcpp::checkUserInterrupt();  // check for ^C from user

        // possible genotypes for this individual
        IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
        int n_poss_gen = poss_gen.size();

        // forward and backward equations
        NumericMatrix alpha = forwardEquations(cross, genotypes(_,ind), founder_geno, is_X_chr, is_female[ind],
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
    if(cross_pu != cross) delete cross_pu;
    delete cross;
    return cur_rec_frac;
}

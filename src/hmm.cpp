// main HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm.h"

// forward equations
NumericMatrix forwardEquations(Cross* cross,
                               IntegerVector genotypes,
                               bool is_X_chr,
                               bool is_female,
                               IntegerVector cross_info,
                               NumericVector rec_frac,
                               IntegerVector marker_index,
                               double error_prob,
                               IntegerVector poss_gen)
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
                                      is_X_chr, is_female, cross_info);
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
                                             is_X_chr, is_female, cross_info);
        }
    }

    return alpha;
}



// backward Equations
NumericMatrix backwardEquations(Cross* cross,
                                IntegerVector genotypes,
                                bool is_X_chr,
                                bool is_female,
                                IntegerVector cross_info,
                                NumericVector rec_frac,
                                IntegerVector marker_index,
                                double error_prob,
                                IntegerVector poss_gen)
{
    int n_pos = marker_index.size();

    // possible genotypes for this chromosome and individual
    int n_gen = poss_gen.size();

    // to contain ln Pr(G_i = g | marker data)
    NumericMatrix beta(n_gen, n_pos);
    
    // initialize beta (not needed; already filled with 0's)
    //    for(int i=0; i<n_gen; i++)
    //        beta(i, n_pos-1) = 0.0;

    // backward equations
    for(int pos = n_pos-2; pos >= 0; pos--) {
        for(int il=0; il<n_gen; il++) {
            beta(il,pos) = beta(0,pos+1) + cross->step(poss_gen[il], poss_gen[0], rec_frac[pos],
                                                       is_X_chr, is_female, cross_info);

            if(marker_index[pos+1] >= 0)
                beta(il,pos) += cross->emit(genotypes[marker_index[pos+1]], poss_gen[0], error_prob, 
                                            is_X_chr, is_female, cross_info);

            for(int ir=1; ir<n_gen; ir++) {
                double to_add = beta(ir,pos+1) + cross->step(poss_gen[il], poss_gen[ir], rec_frac[pos], 
                                                              is_X_chr, is_female, cross_info);
                if(marker_index[pos+1] >=0)
                    to_add += cross->emit(genotypes[marker_index[pos+1]], poss_gen[ir], error_prob,
                                         is_X_chr, is_female, cross_info);
                beta(il,pos) = addlog(beta(il,pos), to_add);
            }
        }
    }

    return beta;
}



// calculate QTL genotype probabilities
// [[Rcpp::export]]
NumericVector calc_genoprob(String crosstype,
                            IntegerMatrix genotypes, // columns are individuals, rows are markers
                            bool is_X_chr,
                            LogicalVector is_female, // length n_ind
                            IntegerMatrix cross_info, // columns are individuals
                            NumericVector rec_frac,   // length nrow(genotypes)-1
                            IntegerVector marker_index, // length nrow(genotypes)
                            double error_prob)
{
    int n_ind = genotypes.cols();
    int n_pos = marker_index.size();

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

    Cross* cross = Cross::Create(crosstype);
    int n_gen = cross->ngen(is_X_chr);
    int matsize = n_gen*n_ind; // size of genotype x individual matrix
    NumericVector genoprobs(matsize*n_pos);

    for(int ind=0; ind<n_ind; ind++) {
        IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
        int n_poss_gen = poss_gen.size();
        NumericMatrix alpha = forwardEquations(cross, genotypes(_,ind), is_X_chr, is_female[ind],
                                               cross_info(_,ind), rec_frac, marker_index, error_prob,
                                               poss_gen);
        NumericMatrix beta = backwardEquations(cross, genotypes(_,ind), is_X_chr, is_female[ind],
                                               cross_info(_,ind), rec_frac, marker_index, error_prob,
                                               poss_gen);

        // calculate genotype probabilities
        for(int pos=0, matindex=n_gen*ind; pos<n_pos; pos++, matindex += matsize) {
            double sum_at_pos = genoprobs[matindex] = alpha(poss_gen[0]-1,pos) + beta(poss_gen[0]-1,pos);
            for(int i=1; i<n_poss_gen; i++) {
                int g = poss_gen[i]-1;
                double val = genoprobs[matindex+g] = alpha(i,pos) + beta(i,pos);
                sum_at_pos = addlog(sum_at_pos, val);
            }
            for(int i=0; i<n_poss_gen; i++) {
                int g = poss_gen[i]-1;
                genoprobs[matindex+g] = exp(genoprobs[matindex+g] - sum_at_pos);
            }
        }
    } // loop over individuals
    genoprobs.attr("dim") = Dimension(n_gen, n_ind, n_pos);
    return genoprobs;
}


/*


// re-estimate inter-marker recombination fractions
double[] estmap(Cross cross,
                GenotypeSymbolMapper[][] genotypes,
                bool is_X_chr,
                bool[] is_female,
                int[][] cross_info,
                Marker[] marker_map,
                double[] rec_frac,
                double error_prob,
                uint max_iterations,
                double tol,
                bool verbose)
{
    if(marker_map.length != rec_frac.length+1)
      throw std::range_error("no. markers in marker map doesn't match rec_frac length");
    if(error_prob < 0.0 || error_prob > 1.0)
      throw std::range_error("error_prob out of range");
    foreach(rf; rec_frac) {
      if(rf < 0 || rf > 0.5)
	throw std::range_error("rec_frac must be >= 0 and <= 0.5");
    }
    if(max_iterations < 0)
      throw std::range_error("max_iterations should be >= 0");
    if(tol < 0)
      throw std::range_error("tol >= 0");
    if(is_female.length != genotypes.length)
      throw std::range_error("is_female should be same length as genotypes");
    if(cross_info.length != genotypes.length)
      throw std::range_error("cross_info should be same length as genotypes");

    auto crossPK = form_cross_phaseknown(cross);

    size_t n_individuals = genotypes.length;
    size_t n_markers = marker_map.length;

    auto cur_rec_frac = rec_frac.dup;
    auto prev_rec_frac = rec_frac.dup;
    auto all_true_geno = crossPK.all_true_geno(is_X_chr);
    double[][] alpha = new double[][](all_true_geno.length, n_markers);
    double[][] beta = new double[][](all_true_geno.length, n_markers);
    double[][] gamma = new double[][](all_true_geno.length, n_markers);
    double sum_gamma;
    foreach(it; 0..max_iterations) {
      foreach(ref rf; cur_rec_frac) {
	rf = 0.0;
      }

      foreach(ind; 0..n_individuals) {

	// forward and backward equations
	alpha = forwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_info[ind], marker_map, prev_rec_frac, error_prob);
	beta = backwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_info[ind], marker_map, prev_rec_frac, error_prob);

        // possible genotypes for this individual, indices to all_true_geno
        auto possible_true_geno_index = crossPK.possible_true_geno_index(is_X_chr, is_female[ind], cross_info[ind]);

	foreach(j; 0..prev_rec_frac.length) {
	  // calculate gamma = log Pr(v1, v2, O)
	  auto sum_gamma_undef = true;
	  foreach(il; possible_true_geno_index) {
            auto left_gen = all_true_geno[il];
	    foreach(ir; possible_true_geno_index) {
              auto right_gen = all_true_geno[ir];
              if(isPseudoMarker(marker_map[j+1]))
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j], is_X_chr, is_female[ind], cross_info[ind]);
              else
                gamma[il][ir] = alpha[il][j] + beta[ir][j+1] +
                  crossPK.emit(genotypes[ind][marker_map[j+1].id], right_gen, error_prob, is_X_chr, is_female[ind], cross_info[ind]) +
                  crossPK.step(left_gen, right_gen, prev_rec_frac[j], is_X_chr, is_female[ind], cross_info[ind]);

	      if(sum_gamma_undef) {
		sum_gamma_undef = false;
		sum_gamma = gamma[il][ir];
	      }
	      else {
		sum_gamma = addlog(sum_gamma, gamma[il][ir]);
	      }
	    }
	  }

	  // update cur_rf
	  foreach(il; possible_true_geno_index) {
	    foreach(ir; possible_true_geno_index) {
	      cur_rec_frac[j] += crossPK.nrec(all_true_geno[il], all_true_geno[ir], is_X_chr, is_female[ind], cross_info[ind]) * exp(gamma[il][ir] - sum_gamma);
	    }
	  }
	} // loop over marker intervals

      } // loop over individuals

      // rescale
      foreach(ref rf; cur_rec_frac) {
	rf /= n_individuals;
	if(rf < tol/1000.0) rf = tol/1000.0;
	else if(rf > 0.5-tol/1000.0) rf = 0.5-tol/1000.0;
      }

      if(verbose) {
	auto maxdif=0.0;
	double tempdif;
	foreach(j; 0..prev_rec_frac.length) {
	  tempdif = abs(prev_rec_frac[j] - cur_rec_frac[j]);
	  if(tempdif > maxdif) {
	    maxdif = tempdif;
	  }
	}
	writefln("%4d %.12f", it, tempdif);
      }

      // check convergence
      auto converged = true;
      foreach(j; 0..prev_rec_frac.length) {
	if(abs(prev_rec_frac[j] - cur_rec_frac[j]) > tol*(cur_rec_frac[j]+tol*100.0)) {
	  converged = false;
	  break;
	}
      }

      if(converged) break;

      prev_rec_frac = cur_rec_frac.dup;
    }

    // calculate log likelihood
    auto loglik = 0.0;
    double curloglik;
    foreach(ind; 0..n_individuals) {

      alpha = forwardEquations(crossPK, genotypes[ind], is_X_chr, is_female[ind], cross_info[ind], marker_map, prev_rec_frac, error_prob);
      auto possible_true_geno_index = crossPK.possible_true_geno_index(is_X_chr, is_female[ind], cross_info[ind]);

      auto curloglik_undef = true;
      foreach(g_index; possible_true_geno_index) {
	if(curloglik_undef) {
	  curloglik_undef = false;
	  curloglik = alpha[g_index][prev_rec_frac.length-1];
	}
	else {
	  curloglik = addlog(curloglik, alpha[g_index][prev_rec_frac.length-1]);
	}
      }
      loglik += curloglik;
    }

    if(verbose) {
      writefln("loglik = %.3f", loglik);
    }

    return(cur_rec_frac);
}

*/

// Calculate addlog(a,b) = log[exp(a) + exp(b)]
// [[Rcpp::export]]
double addlog(const double a, const double b)
{
    const double tol=200.0;

    if(b > a + tol) return(b);
    else if(a > b + tol) return(a);
    else return(a + log1p(exp(b-a)));
}

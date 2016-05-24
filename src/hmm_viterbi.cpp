// Viterbi algoirthm to find arg max Pr(g | O)
//   where g = sequence of true genotypes and O = observed marker genotypes

#include "hmm_viterbi.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "random.h"
#define TOL 1e-6

// find most probable sequence of genotypes
// [[Rcpp::export(".viterbi")]]
IntegerMatrix viterbi(const String& crosstype,
                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                      const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                      const bool is_X_chr,
                      const LogicalVector& is_female, // length n_ind
                      const IntegerMatrix& cross_info, // columns are individuals
                      const NumericVector& rec_frac,   // length nrow(genotypes)-1
                      const IntegerVector& marker_index, // length nrow(genotypes)
                      const double error_prob)
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

    IntegerMatrix result(n_ind, n_pos); // output object

    for(int ind=0; ind<n_ind; ind++) {

        Rcpp::checkUserInterrupt();  // check for ^C from user

        // possible genotypes for this individual
        IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female[ind], cross_info(_,ind));
        int n_poss_gen = poss_gen.size();

        IntegerMatrix traceback(n_pos, n_poss_gen); // for tracing back through the genotypes

        if(n_pos == 1) { // exactly one position
            std::vector<int> best_genotypes;

            // probability of first genotype
            double s = cross->init(poss_gen[0], is_X_chr, is_female[ind], cross_info(_,ind));
            if(marker_index[0] > 0)
                s += cross->emit(genotypes(marker_index[0],ind), poss_gen[0], error_prob,
                                 founder_geno(_, marker_index[0]), is_X_chr, is_female[ind], cross_info(_,ind));
            result(ind,0) = poss_gen[0];

            // probability of other genotypes
            for(int g=1; g<n_poss_gen; g++) {
                double t = cross->init(poss_gen[g], is_X_chr, is_female[ind], cross_info(_,ind));
                if(marker_index[0] > 0)
                    t += cross->emit(genotypes(marker_index[0],ind), poss_gen[g], error_prob,
                                     founder_geno(_, marker_index[0]), is_X_chr, is_female[ind], cross_info(_,ind));
                // bigger or same plus flip coin...bias towards later ones
                if(t > s || (s-t < TOL && R::runif(0.0, 1.0)<0.5)) {
                    s = t;
                    result(ind,0) = poss_gen[g];
                }
            }
        } // exactly one position
        else { // multiple positions
            NumericVector gamma(n_poss_gen);
            NumericVector tempgamma1(n_poss_gen);
            NumericVector tempgamma2(n_poss_gen);

            for(int g=0; g<n_poss_gen; g++) {
                gamma[g] = cross->init(poss_gen[g], is_X_chr, is_female[ind], cross_info(_,ind));
                if(marker_index[0] > 0)
                    gamma[g] += cross->emit(genotypes(marker_index[0],ind), poss_gen[g], error_prob,
                                            founder_geno(_, marker_index[0]), is_X_chr, is_female[ind], cross_info(_,ind));
            }

            for(int pos=0; pos<n_pos-1; pos++) {
                for(int gright=0; gright<n_poss_gen; gright++) {
                    double s = gamma[0] + cross->step(poss_gen[0], poss_gen[gright], rec_frac[pos],
                                                      is_X_chr, is_female[ind], cross_info(_,ind));
                    tempgamma1[gright] = s;
                    traceback(pos,gright) = 0;

                    for(int gleft=1; gleft<n_poss_gen; gleft++) {
                        double t = gamma[gleft] + cross->step(poss_gen[gleft], poss_gen[gright], rec_frac[pos],
                                                              is_X_chr, is_female[ind], cross_info(_,ind));
                        if(t > s || (s-t < TOL && R::runif(0.0, 1.0)<0.5)) {
                            tempgamma1[gright] = s = t;
                            traceback(pos,gright) = gleft;
                        }
                    }
                    if(marker_index[pos+1] >= 0)
                        tempgamma2[gright] = tempgamma1[gright] + cross->emit(genotypes(marker_index[pos+1],ind), poss_gen[gright], error_prob,
                                                                             founder_geno(_, marker_index[pos+1]), is_X_chr, is_female[ind], cross_info(_,ind));
                }
                for(int g=0; g<n_poss_gen; g++) gamma[g] = tempgamma2[g];
            } // loop over positions

            // finish off viterbi and then trace back to get most likely sequence
            result(ind, n_pos-1) = 0;
            double s = gamma[0];
            for(int g=1; g<n_poss_gen; g++) {
                double t = gamma[g];
                if(t > s || (s-t < TOL && R::runif(0.0, 1.0)<0.5)) {
                    s = t;
                    result(ind, n_pos-1) = g;
                }
            }
            for(int pos=n_pos-2; pos>=0; pos--)
                result(ind, pos) = traceback(pos, result(ind, pos+1));

            // replace integers with possible genotypes
            for(int pos=0; pos<n_pos; pos++)
                result(ind,pos) = poss_gen[result(ind,pos)];

        } // if(multiple positions)

    } // loop over individuals

    delete cross;
    return result;
}

// Viterbi algoirthm to find arg max Pr(g | O)
//   where g = sequence of true genotypes and O = observed marker genotypes
// (this version assumes constant is_female and cross_info and pre-calcs the step and emit matrices)

#include "hmm_viterbi2.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "random.h"
#define TOL 1e-6

// find most probable sequence of genotypes
// [[Rcpp::export(".viterbi2")]]
IntegerMatrix viterbi2(const String& crosstype,
                       const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                       const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                       const bool is_X_chr,
                       const bool is_female, // same for all individuals
                       const IntegerVector& cross_info, // same for all individuals
                       const NumericVector& rec_frac,   // length nrow(genotypes)-1
                       const IntegerVector& marker_index, // length nrow(genotypes)
                       const double error_prob)
{
    const int n_ind = genotypes.cols();
    const int n_pos = marker_index.size();
    const int n_mar = genotypes.rows();

    QTLCross* cross = QTLCross::Create(crosstype);

    // check inputs
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
    if(founder_geno.cols() != n_mar)
        throw std::range_error("founder_geno and genotypes have different numbers of markers");
    // end of checks

    IntegerMatrix result(n_ind, n_pos); // output object

    NumericVector init_vector = cross->calc_initvector(is_X_chr, is_female, cross_info);

    int max_obsgeno = max(genotypes);

    std::vector<NumericMatrix> emit_matrix = cross->calc_emitmatrix(error_prob, max_obsgeno, founder_geno,
                                                                    is_X_chr, is_female, cross_info);

    std::vector<NumericMatrix> step_matrix = cross->calc_stepmatrix(rec_frac, is_X_chr, is_female, cross_info);

    // possible genotypes
    IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female, cross_info);
    int n_poss_gen = poss_gen.size();

    for(int ind=0; ind<n_ind; ind++) {

        Rcpp::checkUserInterrupt();  // check for ^C from user

        IntegerMatrix traceback(n_pos, n_poss_gen); // for tracing back through the genotypes

        if(n_pos == 1) { // exactly one position
            std::vector<int> best_genotypes;

            // probability of first genotype
            double s = init_vector[0];
            if(marker_index[0] >= 0)
                s += emit_matrix[marker_index[0]](genotypes(marker_index[0],ind), 0);
            result(ind,0) = poss_gen[0];

            // probability of other genotypes
            for(int g=1; g<n_poss_gen; g++) {
                double t = init_vector[g];
                if(marker_index[0] >= 0)
                    t += emit_matrix[marker_index[0]](genotypes(marker_index[0],ind), g);
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
                gamma[g] = init_vector[g];
                if(marker_index[0] >= 0)
                    gamma[g] += emit_matrix[marker_index[0]](genotypes(marker_index[0],ind), g);
            }

            for(int pos=0; pos<n_pos-1; pos++) {
                for(int gright=0; gright<n_poss_gen; gright++) {
                    double s = gamma[0] + step_matrix[pos](0, gright);
                    tempgamma1[gright] = s;
                    traceback(pos,gright) = 0;

                    for(int gleft=1; gleft<n_poss_gen; gleft++) {
                        double t = gamma[gleft] + step_matrix[pos](gleft, gright);
                        if(t > s || (s-t < TOL && R::runif(0.0, 1.0)<0.5)) {
                            tempgamma1[gright] = s = t;
                            traceback(pos,gright) = gleft;
                        }
                    }
                    if(marker_index[pos+1] >= 0)
                        tempgamma2[gright] = tempgamma1[gright] + emit_matrix[marker_index[pos+1]](genotypes(marker_index[pos+1],ind), gright);
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

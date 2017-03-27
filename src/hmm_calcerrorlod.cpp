// calculate genotyping error LOD scores
// (assumes constant is_female and cross_info and pre-calcs the step and emit matrices)

#include "hmm_calcerrorlod.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "hmm_util.h"
#include "hmm_forwback2.h"

// calculate genotyping error lod scores (output is mar x ind and so should be transposed)
// [[Rcpp::export(".calc_errorlod")]]
NumericMatrix calc_errorlod(const String& crosstype,
                            const NumericVector& probs, // genotype probs [genotype, ind, marker]
                            const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                            const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                            const bool is_X_chr,
                            const bool is_female, // same for all individuals
                            const IntegerVector& cross_info) // same for all individuals
{
    const double error_prob = 0.01; // just used to get emit values, to determine errors from non-errors

    const int n_ind = genotypes.cols();
    const int n_mar = genotypes.rows();
    if(Rf_isNull(probs.attr("dim")))
        throw std::invalid_argument("probs has no dimension attribute");
    const IntegerVector& dim_probs = probs.attr("dim");
    if(dim_probs.size() != 3)
        throw std::invalid_argument("probs should be 3-dimensional array of probabilities");
    const double log_half = log(0.5);

    QTLCross* cross = QTLCross::Create(crosstype);

    // check inputs
    if(n_ind != dim_probs[1])
        throw std::invalid_argument("Different number of individuals in genotypes and probabilities.");
    if(n_mar != dim_probs[2])
        throw std::invalid_argument("Different number of markers in genotypes and probabilities.");
    if(!cross->check_founder_geno_size(founder_geno, n_mar))
        throw std::range_error("founder_geno is not the right size");
    if(founder_geno.cols() != n_mar)
        throw std::range_error("founder_geno and genotypes have different numbers of markers");
    // end of checks

    const int n_gen = cross->ngen(is_X_chr);
    NumericMatrix error_lod(n_mar, n_ind);

    NumericVector init_vector = cross->calc_initvector(is_X_chr, is_female, cross_info);

    const int matsize = n_ind * n_gen;
    const int max_obsgeno = max(genotypes);

    std::vector<NumericMatrix> emit_matrix = cross->calc_emitmatrix(error_prob, max_obsgeno,
                                                                    founder_geno,
                                                                    is_X_chr, is_female, cross_info);

    // possible genotypes
    IntegerVector poss_gen = cross->possible_gen(is_X_chr, is_female, cross_info);
    const int n_poss_gen = poss_gen.size();

    for(int ind=0; ind<n_ind; ind++) {

        Rcpp::checkUserInterrupt();  // check for ^C from user

        for(int mar=0, matindex=n_gen*ind; mar<n_mar; mar++, matindex += matsize) {

            double init_err=0.0, init_noerr=0.0, post_err=0.0, post_noerr=0.0;

            int obs_geno = genotypes(mar,ind);
            if(obs_geno == 0) { // missing genotype
                error_lod(mar, ind) = 0.0;
                continue;
            }

            for(int i=0; i<n_poss_gen; i++) {
                int g = poss_gen[i]-1;
                if(emit_matrix[mar](obs_geno, i) < log_half) { // error
                    init_err += exp(init_vector[i]);
                    post_err += probs[matindex + g];
                }
                else { // not an error
                    init_noerr += exp(init_vector[i]);
                    post_noerr += probs[matindex + g];
                }
            } // loop over possible genotypes

            error_lod(mar, ind) = log10( (post_err * init_noerr) / (post_noerr * init_err) );

        } // loop over markers
    } // loop over individuals

    return error_lod;
}

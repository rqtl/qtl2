// phase-known AIL3 (3-way advanced intercross lines) QTLCross class (for HMM, in particular est.map)
// (assuming all F1 hybrids formed followed by random mating with large population)

#ifndef CROSS_AIL3PK_H
#define CROSS_AIL3PK_H

#include <Rcpp.h>
#include "cross.h"

class AIL3PK : public QTLCross
{
 public:
    AIL3PK(){
        crosstype = "ail3pk";
        phase_known_crosstype = "ail3pk";
    };
    ~AIL3PK(){};

    const bool check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const double init(const int true_gen,
                      const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);
    const double emit(const int obs_gen, const int true_gen, const double error_prob,
                      const Rcpp::IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const Rcpp::IntegerVector& cross_info);
    const double step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const Rcpp::IntegerVector possible_gen(const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    const int ngen(const bool is_x_chr);

    const int nalleles();

    const int nrec(const int gen_left, const int gen_right,
                   const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

    // this isn't actually implemented yet; just throws an error
    const double est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                              const Rcpp::IntegerMatrix& cross_info, const int n_gen);

    // this is to indicate that the f2pk cross type shouldn't exist on the R side
    // (it's strictly a device for using phase-known version for est_map)
    const bool crosstype_supported() {
        return false;
    }

    const Rcpp::NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const Rcpp::LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr);

};

#endif // CROSS_AIL3PK_H

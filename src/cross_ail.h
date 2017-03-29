// AIL QTLCross class (for HMM)

#ifndef CROSS_AIL_H
#define CROSS_AIL_H

#include <Rcpp.h>
#include "cross.h"

class AIL : public QTLCross
{
 public:
    AIL(){
        crosstype = "ail";
        phase_known_crosstype = "ailpk";
    };
    ~AIL(){};

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

    const Rcpp::NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const Rcpp::LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr);

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles, const bool is_x_chr);

    const int nrec(const int gen_left, const int gen_right,
                   const bool is_x_chr, const bool is_female,
                   const Rcpp::IntegerVector& cross_info);
};

#endif // CROSS_AIL_H

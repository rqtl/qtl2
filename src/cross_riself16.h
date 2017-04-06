// 16-way RIL by selfing QTLCross class (for HMM)

#ifndef CROSS_RISELF16_H
#define CROSS_RISELF16_H

#include <Rcpp.h>
#include "cross.h"

class RISELF16 : public QTLCross
{
 public:
    RISELF16(){
        crosstype = "riself16";
        phase_known_crosstype = "riself16";
    };
    ~RISELF16(){};

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

    const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr);

    const bool check_founder_geno_size(const Rcpp::IntegerMatrix& founder_geno, const int n_markers);
    const bool check_founder_geno_values(const Rcpp::IntegerMatrix& founder_geno);
    const bool need_founder_geno();

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles, const bool is_x_chr);

    const int nrec(const int gen_left, const int gen_right,
                   const bool is_x_chr, const bool is_female,
                   const Rcpp::IntegerVector& cross_info);

    const double est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                              const Rcpp::IntegerMatrix& cross_info, const int n_gen);

    // check whether X chr can be handled
    const bool check_handle_x_chr(const bool any_x_chr);

};

#endif // CROSS_RISELF16_H

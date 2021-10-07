// Heterogeneous stock F1 (in cross with another inbred strain) QTLCross class (for HMM)

#ifndef CROSS_HSF1_H
#define CROSS_HSF1_H

#include <Rcpp.h>
#include "cross.h"

class HSF1 : public QTLCross
{
 public:
    HSF1(){
        crosstype = "hsf1";
        phase_known_crosstype = "hsf1";
    };
    ~HSF1(){};

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

    const Rcpp::NumericMatrix geno2allele_matrix(const bool is_x_chr);

    const bool check_is_female_vector(const Rcpp::LogicalVector& is_female, const bool any_x_chr);

    const bool check_crossinfo(const Rcpp::IntegerMatrix& cross_info, const bool any_x_chr);

    const bool check_founder_geno_size(const Rcpp::IntegerMatrix& founder_geno, const int n_markers);
    const bool check_founder_geno_values(const Rcpp::IntegerMatrix& founder_geno);
    const bool need_founder_geno();

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles, const bool is_x_chr);

    const double est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                              const Rcpp::IntegerMatrix& cross_info, const int n_gen);

    // not yet implemented
    const Rcpp::NumericVector est_map2(const Rcpp::IntegerMatrix& genotypes,
                                       const Rcpp::IntegerMatrix& founder_geno,
                                       const bool is_X_chr,
                                       const Rcpp::LogicalVector& is_female,
                                       const Rcpp::IntegerMatrix& cross_info,
                                       const Rcpp::IntegerVector& cross_group,
                                       const Rcpp::IntegerVector& unique_cross_group,
                                       const Rcpp::NumericVector& rec_frac,
                                       const double error_prob,
                                       const int max_iterations,
                                       const double tol,
                                       const bool verbose);

};

#endif // CROSS_HSF1_H

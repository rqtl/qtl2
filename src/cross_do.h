// Diversity Outcross QTLCross class (for HMM)

#ifndef CROSS_DO_H
#define CROSS_DO_H

#include <Rcpp.h>
#include "cross.h"

class DO : public QTLCross
{
 public:
    DO(){
        crosstype = "do";
        phase_known_crosstype = "dopk";
    };
    ~DO(){};

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

    const bool is_het(const int true_gen); // is heterozygous (just for female X or autosome)

    const int encode_alleles(const int allele1, const int allele2); // convert (a1,a2) pair to genotype 1-36
    const Rcpp::IntegerVector decode_geno(const int true_gen);      // convert genotype to (a1,a2) pair

    // helper functions for step()
    const double step_auto(int left, int right, double r, int s,
                           Rcpp::IntegerVector precc_gen, Rcpp::NumericVector precc_alpha);
    const double step_femX(int left, int right, double r, int s,
                           Rcpp::IntegerVector precc_gen, Rcpp::NumericVector precc_alpha);
    const double step_malX(int left, int right, double r, int s,
                           Rcpp::IntegerVector precc_gen, Rcpp::NumericVector precc_alpha);

    const bool check_founder_geno_size(const Rcpp::IntegerMatrix& founder_geno, const int n_markers);
    const bool check_founder_geno_values(const Rcpp::IntegerMatrix& founder_geno);
    const bool need_founder_geno();

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles, const bool is_x_chr);
};

#endif // CROSS_DO_H

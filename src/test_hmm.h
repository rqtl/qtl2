// functions to test basic HMM things from R
#ifndef TEST_HMM_H
#define TEST_HMM_H

#include <Rcpp.h>

// test init functions from R
double test_init(const Rcpp::String& crosstype,
                 const int true_gen,
                 const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

// test emit functions from R
double test_emit(const Rcpp::String& crosstype,
                 const int obs_gen, const int true_gen, const double error_prob,
                 const Rcpp::IntegerVector& founder_geno, const bool is_x_chr,
                 const bool is_female, const Rcpp::IntegerVector& cross_info);

// test emit functions from R
double test_step(const Rcpp::String& crosstype,
                 const int gen_left, const int gen_right, const double rec_frac,
                 const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

bool test_check_geno(const Rcpp::String& crosstype, const int gen, const bool is_observed_value,
                     const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

Rcpp::IntegerVector test_possible_gen(const Rcpp::String& crosstype,
                                      const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

int test_ngen(const Rcpp::String& crosstype, const bool is_x_chr);

double test_nrec(const Rcpp::String& crosstype, const int gen_left, const int gen_right,
                 const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

bool test_founder_geno(const Rcpp::String& crosstype, const Rcpp::IntegerMatrix& founder_geno, const int n_markers);

bool need_founder_geno(const Rcpp::String& crosstype);

// test calculation of vector of emit matrices
std::vector<Rcpp::NumericMatrix> test_emitmatrix(const Rcpp::String& crosstype,
                                                 const double error_prob,
                                                 const int max_obsgeno,
                                                 const Rcpp::IntegerMatrix& founder_geno, const bool is_x_chr,
                                                 const bool is_female, const Rcpp::IntegerVector& cross_info);

// test calculation of vector of transition matrices
std::vector<Rcpp::NumericMatrix> test_stepmatrix(const Rcpp::String& crosstype,
                                                 const Rcpp::NumericVector& rec_frac,
                                                 const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

// test calculation of init vector
Rcpp::NumericVector test_initvector(const Rcpp::String& crosstype,
                                    const bool is_x_chr, const bool is_female, const Rcpp::IntegerVector& cross_info);

#endif // TEST_HMM_H

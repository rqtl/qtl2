// random number generation (e.g., permutations)
#ifndef RANDOM_H
#define RANDOM_H

#include <vector>
#include <map>
#include <Rcpp.h>

// random integer from {low, low+1, ..., high}
int random_int(const int low, const int high);

// vector of random integers from {low, low+1, ..., high}
Rcpp::IntegerVector random_int(const int n, const int low, const int high);

// permute a vector of numbers
Rcpp::NumericVector permute_nvector(const Rcpp::NumericVector x);
Rcpp::IntegerVector permute_ivector(const Rcpp::IntegerVector x);
std::vector<double> permute_nvector(const std::vector<double> x);
std::vector<int> permute_ivector(const std::vector<int> x);

// permute a vector of numbers in place
void permute_nvector_inplace(Rcpp::NumericVector x);
void permute_ivector_inplace(Rcpp::IntegerVector x);
void permute_nvector_inplace(std::vector<double> x);
void permute_ivector_inplace(std::vector<int> x);

// get permutation of {0..(n-1)}
Rcpp::IntegerVector get_permutation(const int n);

// get a set of permutations of a vector, as columns of a matrix
Rcpp::NumericMatrix permute_nvector(const int n_perm,
                                    const Rcpp::NumericVector x);
Rcpp::IntegerMatrix permute_ivector(const int n_perm,
                                    const Rcpp::IntegerVector x);

// stratified permutation
Rcpp::NumericMatrix permute_nvector_stratified(const int n_perm,
                                               const Rcpp::NumericVector& x,
                                               const Rcpp::IntegerVector& strata,
                                               int n_strata);
Rcpp::IntegerMatrix permute_ivector_stratified(const int n_perm,
                                               const Rcpp::IntegerVector& x,
                                               const Rcpp::IntegerVector& strata,
                                               int n_strata);

#endif // RANDOM_H

// calculate matrix of counts of genotype matches for pairs of individuals
#ifndef COMPARE_GENO_H
#define COMPARE_GENO_H

#include <Rcpp.h>

Rcpp::IntegerMatrix compare_geno(const Rcpp::IntegerMatrix& geno); // matrix n_mar x n_ind (transposed of normal)

#endif // COMPARE_GENO_H

// guess phase in imputed genotypes
#ifndef GUESS_PHASE_H
#define GUESS_PHASE_H


#include <Rcpp.h>

// guess phase, F2 autosome
Rcpp::IntegerVector guess_phase_f2A(const Rcpp::IntegerMatrix& geno,
                                    bool deterministic); // pos x ind


// guess phase, F2 X chr
Rcpp::IntegerVector guess_phase_f2X(const Rcpp::IntegerMatrix& geno,
                                    bool deterministic); // pos x ind

// guess phase, MPP autosome
Rcpp::IntegerVector guess_phase_A(const Rcpp::IntegerMatrix& geno, // pos x ind
                                  const Rcpp::String& crosstype,
                                  bool deterministic);

// guess phase, MPP X chr
Rcpp::IntegerVector guess_phase_X(const Rcpp::IntegerMatrix& geno, // pos x ind
                                  const Rcpp::String& crosstype,
                                  const Rcpp::LogicalVector& is_female,
                                  bool deterministic);


// guess phase of genotypes for one individual along one chromosome
Rcpp::IntegerVector phase_geno(Rcpp::IntegerVector g1, Rcpp::IntegerVector g2,
                               bool deterministic);

#endif // GUESS_PHASE_H

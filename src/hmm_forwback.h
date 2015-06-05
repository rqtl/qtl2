// forward-backward equations for HMM
#ifndef HMM_FORWBACK_H
#define HMM_FORWBACK_H

#include "cross.h"
#include "hmm_util.h"

// forward equations
NumericMatrix forwardEquations(QTLCross* cross,
                               const IntegerVector& genotypes,
                               const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                               const bool is_X_chr,
                               const bool is_female,
                               const IntegerVector& cross_info,
                               const NumericVector& rec_frac,
                               const IntegerVector& marker_index,
                               const double error_prob,
                               const IntegerVector& poss_gen);


// backward Equations
NumericMatrix backwardEquations(QTLCross* cross,
                                const IntegerVector& genotypes,
                                const IntegerMatrix& founder_geno, // columns are markers, rows are founder lines
                                const bool is_X_chr,
                                const bool is_female,
                                const IntegerVector& cross_info,
                                const NumericVector& rec_frac,
                                const IntegerVector& marker_index,
                                const double error_prob,
                                const IntegerVector& poss_gen);

#endif // HMM_FORWBACK_H

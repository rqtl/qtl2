// general qtlcross class
//
// see cross.cpp for info on how to add a new cross type

#ifndef CROSS_H
#define CROSS_H

#include <Rcpp.h>

using namespace Rcpp;

class QTLCross
{
public:
    String type;

    String phase_known_type;

    static QTLCross* Create(String type);

    virtual bool check_geno(int gen, bool is_observed_value,
                            bool is_X_chr, bool is_female,
                            IntegerVector cross_info) {
        return false;
    }

    virtual double init(int true_gen,
                        bool is_X_chr, bool is_female,
                        IntegerVector cross_info) {
        return 0.0;
    }

    virtual double emit(int obs_gen, int true_gen, double error_prob,
                        bool is_X_chr, bool is_female,
                        IntegerVector cross_info) {
        return 0.0;
    }

    virtual double step(int gen_left, int gen_right, double rec_frac,
                        bool is_X_chr, bool is_female,
                        IntegerVector cross_info) {
        return 0.0;
    }

    virtual int ngen(bool is_X_chr) {
        return 0;
    }

    virtual IntegerVector possible_gen(bool is_X_chr, bool is_female,
                                       IntegerVector cross_info) {
        int ng = ngen(is_X_chr);
        IntegerVector x(ng);
        for(int i=0; i<ng; i++) x[i] = i+1;
        return x;
    }

    virtual double nrec(int gen_left, int gen_right,
                        bool is_X_chr, bool is_female,
                        IntegerVector cross_info) {
        return 0.0;
    }

    virtual double est_rec_frac(NumericMatrix full_gamma, bool is_X_chr) {
        return 0.0;
    }

};

#endif // CROSS.H

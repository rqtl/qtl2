// general cross class + cross factory
//
// to add a new cross type:
//     - create files similar to cross_f2.h and cross_f2.cpp
//     - add include line below
//     - add if statement within Cross::Create function below
//
// to create a Cross instance using a string with cross type:
//     Cross* cross = Cross::Create("f2");
// then refer to functions like cross->init()

#ifndef CROSS_H
#define CROSS_H

#include <Rcpp.h>

using namespace Rcpp;

class Cross
{
public:
    String type;

    String phase_known_type;

    static Cross* Create(String type);

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

};

#include "cross_f2.h"
#include "cross_f2pk.h"
#include "cross_bc.h"
#include "cross_risib.h"
#include "cross_riself.h"

#endif // CROSS.H

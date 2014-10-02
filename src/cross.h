// general cross class + cross factory
//
// to add a new cross type:
//     - create files similar to cross_f2.h and cross_f2.cpp
//     - add include line below
//     - add if statement within Cross::Create function below
//
// to create a Cross instance using a string with cross type:
//     Cross* cross1 = Cross::Create("f2");
// or  Cross  cross1 = *(Cross::Create("f2"));

#ifndef CROSS_H
#define CROSS_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Cross
{
public:
    string type;

    static Cross* Create(string type);

    virtual double init(int true_gen,
                        bool is_X_chr, bool is_female,
                        vector<int>cross_info) {
        return 0.0;
    }

    virtual double emit(int obs_gen, int true_gen, double error_prob,
                        bool is_X_chr, bool is_female,
                        vector<int> cross_info) {
        return 0.0;
    }

    virtual double step(int gen_left, int gen_right, double rec_frac,
                        bool is_X_chr, bool is_female,
                        vector<int> cross_info) {
        return 0.0;
    }

    virtual vector<int> geno(bool is_X_chr, bool is_female,
                             vector<int> cross_info) {
        vector<int> x(0);
        return x;
    }

    virtual double nrec(int gen_left, int gen_right,
                        bool is_X_chr, bool is_female,
                        vector<int> cross_info) {
        return 0.0;
    }

    virtual double initPK(int true_gen,
                        bool is_X_chr, bool is_female,
                        vector<int>cross_info) {
        return init(true_gen, is_X_chr, is_female, cross_info);
    }

    virtual double emitPK(int obs_gen, int true_gen, double error_prob,
                        bool is_X_chr, bool is_female,
                        vector<int> cross_info) {
        return emit(obs_gen, true_gen, error_prob, is_X_chr, is_female, cross_info);
    }

    virtual double stepPK(int gen_left, int gen_right, double rec_frac,
                        bool is_X_chr, bool is_female,
                        vector<int> cross_info) {
        return step(gen_left, gen_right, rec_frac, is_X_chr, is_female, cross_info);
    }

    virtual vector<int> genoPK(bool is_X_chr, bool is_female,
                               vector<int> cross_info) {
        return geno(is_X_chr, is_female, cross_info);
    }
};

#include "cross_f2.h"
#include "cross_bc.h"
#include "cross_risib.h"
#include "cross_riself.h"

#endif // CROSS.H

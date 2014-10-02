#ifndef CROSS_RISELF_H
#define CROSS_RISELF_H

#include <vector>
#include <math.h>
#include "cross.h"

class RIself : public Cross
{
 public:
    RIself(){
        type = "riself";
    };
    ~RIself(){};

    double init(int true_gen, bool is_X_chr, bool is_female,
                vector<int> cross_info);
    double emit(int obs_gen, int true_gen, double error_prob, 
                bool is_X_chr, bool is_female, vector<int> cross_info);
    double step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female, vector<int> cross_info);
    vector<int> geno(bool is_X_chr, bool is_female,
                     vector <int>cross_info);
};

#endif // CROSS_RISELF_H

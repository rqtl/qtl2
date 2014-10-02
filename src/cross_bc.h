#ifndef CROSS_BC_H
#define CROSS_BC_H

#include <vector>
#include <math.h>
#include "cross.h"

class BC : public Cross
{
 public:
    BC(){
        type = "bc";
    };
    ~BC(){};

    bool check_geno(int gen, bool is_observed_value,
                    bool is_X_chr, bool is_female, vector<int> cross_info);

    double init(int true_gen, bool is_X_chr, bool is_female,
                vector<int> cross_info);
    double emit(int obs_gen, int true_gen, double error_prob, 
                bool is_X_chr, bool is_female, vector<int> cross_info);
    double step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female, vector<int> cross_info);

    vector<int> geno(bool is_X_chr, bool is_female,
                     vector <int>cross_info);
    vector<int> allgeno(bool is_X_chr);

    double nrec(int gen_left, int gen_right,
                bool is_X_chr, bool is_female, vector<int> cross_info);
};

#endif // CROSS_BC_H

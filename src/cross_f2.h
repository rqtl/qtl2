#ifndef CROSS_F2_H
#define CROSS_F2_H

#include <vector>
#include <math.h>
#include "cross.h"

using namespace std;

class F2 : public Cross
{
 public:
    F2(){
        type = "f2";
    };
    ~F2(){};

    double init(int true_gen, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info);
    double emit(int obs_gen, int true_gen, double error_prob, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info);
    double step(int gen_left, int gen_right, double rec_frac, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info);
};

#endif // CROSS_F2_H

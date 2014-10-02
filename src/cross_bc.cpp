// backcross HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"

enum ogen {NA=0, A=1, B=2};
enum tgen {AA=1, AB=2, AY=1, BY=3};

double BC::init(int true_gen, bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    return log(0.5);
}

double BC::emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(obs_gen) {
    case NA: return 0.0; // missing value
    case AA: case AB:
        if(obs_gen==true_gen) return log(1.0-error_prob);
        else return log(error_prob);
    }

    Rcpp::exception("invalid obs_gen or true_gen");
    return NA_REAL; // can't get here
}


double BC::step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    if(gen_left==gen_right) return log(1.0-rec_frac);
    else return log(rec_frac);

    Rcpp::exception("invalid genotypes");
    return NA_REAL; // can't get here
}

vector<int> BC::geno(bool is_X_chr, bool is_female,
                 vector <int>cross_info)
{
    int vals[] = {AA,AB};
    vector<int> result(vals, vals+2);
    return result;
}

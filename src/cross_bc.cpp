#include <math.h>
#include <Rcpp.h>
#include "cross.h"

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
    case 0: return(0.0); // missing value
    case 1: case 2:
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
    int vals[] = {1,2};
    vector<int> result(vals, vals+2);
    return result;
}

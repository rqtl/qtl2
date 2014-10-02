#include <math.h>
#include <Rcpp.h>
#include "cross.h"

double F2::init(int true_gen,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    if(true_gen==2) return log(0.5);
    if(true_gen==1 || true_gen==3) return log(0.25);

    Rcpp::exception("invalid genotype");
    return NA_REAL; // can't get here
}

double F2::emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(obs_gen) {
    case 0: return(0.0);
    case 1: case 2: case 3:
        if(obs_gen==true_gen) return log(1.0-error_prob);
        else return log(error_prob*0.5);
    case 4: /* AA or AB (not BB) */
        if(true_gen != 3) return log(1.0-error_prob*0.5);
        else return log(error_prob);
    case 5: /* AB or BB (not AA) */
        if(true_gen != 1) return log(1.0-error_prob*0.5);
        else return log(error_prob);
    }
    Rcpp::exception("invalid obs_gen or true_gen");
    return NA_REAL; /* shouldn't get here */
}


double F2::step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(gen_left) {
    case 1:
        switch(gen_right) {
        case 1: return 2.0*log(1.0-rec_frac);
        case 2: return log(0.5)+log(1.0-rec_frac)+log(rec_frac);
        case 3: return 2.0*log(rec_frac);
        }
    case 2:
        switch(gen_right) {
        case 1: case 3: return log(rec_frac)+log(1.0-rec_frac);
        case 2: return log((1.0-rec_frac)*(1.0-rec_frac)+rec_frac*rec_frac);
        }
    case 3:
        switch(gen_right) {
        case 1: return 2.0*log(rec_frac);
        case 2: return log(0.5)+log(1.0-rec_frac)+log(rec_frac);
        case 3: return 2.0*log(1.0-rec_frac);
        }
    }
    Rcpp::exception("invalid genotypes.");
    return NA_REAL; // can't get here
}

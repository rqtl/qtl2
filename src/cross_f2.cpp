// intercross HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"

enum ogen {NA=0, A=1, H=2, B=3, notB=4, notA=5};
enum tgen {AA=1, AB=2, BB=3, AY=1, BY=3};

double F2::init(int true_gen,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    if(true_gen==AB) return log(0.5);
    if(true_gen==AA || true_gen==BB) return log(0.25);

    Rcpp::exception("invalid genotype");
    return NA_REAL; // can't get here
}

double F2::emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(obs_gen) {
    case NA: return 0.0;
    case A: case H: case B:
        if(obs_gen==true_gen) return log(1.0-error_prob);
        else return log(error_prob*0.5);
    case notB: /* AA or AB (not BB) */
        if(true_gen != BB) return log(1.0-error_prob*0.5);
        else return log(error_prob);
    case notA: /* AB or BB (not AA) */
        if(true_gen != AA) return log(1.0-error_prob*0.5);
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
    case AA:
        switch(gen_right) {
        case AA: return 2.0*log(1.0-rec_frac);
        case AB: return log(0.5)+log(1.0-rec_frac)+log(rec_frac);
        case BB: return 2.0*log(rec_frac);
        }
    case AB:
        switch(gen_right) {
        case AA: case BB: return log(rec_frac)+log(1.0-rec_frac);
        case AB: return log((1.0-rec_frac)*(1.0-rec_frac)+rec_frac*rec_frac);
        }
    case BB:
        switch(gen_right) {
        case AA: return 2.0*log(rec_frac);
        case AB: return log(0.5)+log(1.0-rec_frac)+log(rec_frac);
        case BB: return 2.0*log(1.0-rec_frac);
        }
    }

    Rcpp::exception("invalid genotypes.");
    return NA_REAL; // can't get here
}

vector<int> F2::geno(bool is_X_chr, bool is_female,
                 vector <int>cross_info)
{
    int vals[] = {AA,AB,BB};
    vector<int> result(vals, vals+3);
    return result;
}

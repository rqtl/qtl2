#include <math.h>
#include <string>
#include <Rcpp.h>
#include "cross.h"

double RIself::init(int true_gen, bool ignored1, bool ignored2,
                    vector<int> ignored3)
{
    if(true_gen == 1 || true_gen == 2) return -log(2.0);

    Rcpp::exception("invalid genotype");
    return NA_REAL;
}

double RIself::emit(int obs_gen, int true_gen, double error_prob,
                    bool ignored1, bool ignored2,
                    vector<int> ignored3)
{
    if(obs_gen==0) return 0.0; // missing

    if(obs_gen!=1 && obs_gen!=2)
        Rcpp::exception("Incorrect value for obs_gen");
    if(true_gen!=1 && true_gen!=2)
        Rcpp::exception("Incorrect value for true_gen");

    if(obs_gen == true_gen) return log(1.0 - error_prob);
    else return log(error_prob);
}

double RIself::step(int gen_left, int gen_right, double rec_frac,
                    bool ignored1, bool ignored2,
                    vector<int> ignored3)
{
    const double R = 2.0*rec_frac/(1+2.0*rec_frac);

    if(gen_left!=1 && gen_left!=2)
        Rcpp::exception("Incorrect value for gen_left");
    if(gen_right!=1 && gen_right!=2)
        Rcpp::exception("Incorrect value for gen_right");

    if(gen_left == gen_right) return log(1.0-R);
    else return log(R);
}

vector<int> RIself::geno(bool ignored1, bool ignored2,
                         vector<int> ignored3)
{
    int vals[] = {1,2};
    vector<int> result(vals, vals+2);
    return result;
}

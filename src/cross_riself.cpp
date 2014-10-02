// RI selfing HMM functions

#include <math.h>
#include <string>
#include <Rcpp.h>
#include "cross.h"

enum ogen {NA=0, A=1, B=2};
enum tgen {AA=1, BB=2};

double RIself::init(int true_gen, bool ignored1, bool ignored2,
                    vector<int> ignored3)
{
    if(true_gen == AA || true_gen == BB) return -log(2.0);

    Rcpp::exception("invalid genotype");
    return NA_REAL;
}

double RIself::emit(int obs_gen, int true_gen, double error_prob,
                    bool ignored1, bool ignored2, vector<int> ignored3)
{
    if(obs_gen==NA) return 0.0; // missing

    if(obs_gen!=A && obs_gen!=B)
        Rcpp::exception("Incorrect value for obs_gen");
    if(true_gen!=AA && true_gen!=BB)
        Rcpp::exception("Incorrect value for true_gen");

    if(obs_gen == true_gen) return log(1.0 - error_prob);
    else return log(error_prob);
}

double RIself::step(int gen_left, int gen_right, double rec_frac,
                    bool ignored1, bool ignored2, vector<int> ignored3)
{
    const double R = 2.0*rec_frac/(1+2.0*rec_frac);

    if(gen_left!=AA && gen_left!=BB)
        Rcpp::exception("Incorrect value for gen_left");
    if(gen_right!=AA && gen_right!=BB)
        Rcpp::exception("Incorrect value for gen_right");

    if(gen_left == gen_right) return log(1.0-R);
    else return log(R);
}

vector<int> RIself::geno(bool ignored1, bool ignored2, vector<int> ignored3)
{
    int vals[] = {AA,BB};
    vector<int> result(vals, vals+2);
    return result;
}

double RIself::nrec(int gen_left, int gen_right,
                    bool ignored1, bool ignored2, vector<int> ignored3)
{
    if(gen_left!=AA && gen_left!=BB)
        Rcpp::exception("Incorrect value for gen_left");
    if(gen_right!=AA && gen_right!=BB)
        Rcpp::exception("Incorrect value for gen_right");

    if(gen_left == gen_right) return 0.0;
    else return 1.0;
}

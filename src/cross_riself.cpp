// RI selfing HMM functions

#include "cross.h"

enum gen {GENO_NA=0, AA=1, BB=2};

bool RIself::check_geno(int gen, bool is_observed_value,
                        bool ignored1, bool ignored2, IntegerVector ignored3)
{
    if(is_observed_value && gen==GENO_NA) return true;

    if(gen==AA || gen==BB) return(true);

    throw std::range_error("invalid genotype");
    return false; // can't get here
}

double RIself::init(int true_gen, bool ignored1, bool ignored2,
                    IntegerVector ignored3)
{
    check_geno(true_gen, false, ignored1, ignored2, ignored3);

    return -log(2.0);
}

double RIself::emit(int obs_gen, int true_gen, double error_prob,
                    bool ignored1, bool ignored2, IntegerVector ignored3)
{
    check_geno(obs_gen, true, ignored1, ignored2, ignored3);
    check_geno(true_gen, false, ignored1, ignored2, ignored3);

    if(obs_gen==GENO_NA) return 0.0; // missing

    if(obs_gen == true_gen) return log(1.0 - error_prob);
    else return log(error_prob);
}

double RIself::step(int gen_left, int gen_right, double rec_frac,
                    bool ignored1, bool ignored2, IntegerVector ignored3)
{
    check_geno(gen_left, false, ignored1, ignored2, ignored3);
    check_geno(gen_right, false, ignored1, ignored2, ignored3);

    const double R = 2.0*rec_frac/(1+2.0*rec_frac);

    if(gen_left == gen_right) return log(1.0-R);
    else return log(R);
}

IntegerVector RIself::allgeno(bool ignored1)
{
    int vals[] = {AA,BB};
    IntegerVector result(vals, vals+2);
    return result;
}

double RIself::nrec(int gen_left, int gen_right,
                    bool ignored1, bool ignored2, IntegerVector ignored3)
{
    check_geno(gen_left, false, ignored1, ignored2, ignored3);
    check_geno(gen_right, false, ignored1, ignored2, ignored3);

    if(gen_left == gen_right) return 0.0;
    else return 1.0;
}

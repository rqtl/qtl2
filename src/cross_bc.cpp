// backcross HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_bc.h"

enum gen {AA=1, AB=2, BB=3, AY=3, BY=4};

bool BC::check_geno(int gen, bool is_observed_value,
                    bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    if(is_observed_value && gen==0) return true;

    if(!is_X_chr || (is_X_chr && is_female)) {
        if(gen != AA && gen != AB)
            throw std::range_error("genotype value not allowed");
    }
    else { // male X chr
        if(is_observed_value) {
            if(gen != AA && gen != BB)
                throw std::range_error("genotype value not allowed");
        }
        else {
            if(gen != AY && gen != BY)
                throw std::range_error("genotype value not allowed");
        }
    }

    return true;
}

double BC::init(int true_gen,
                bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    return log(0.5);
}

double BC::emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    check_geno(obs_gen, true, is_X_chr, is_female, cross_info);
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    if(obs_gen==0) return 0.0; // missing value

    if(is_X_chr && !is_female) { // X chr males different
        if(obs_gen==AA) {
            if(true_gen==AY) return log(1.0-error_prob);
            else return log(error_prob);
        }
        else { // BB
            if(true_gen==BY) return log(1.0-error_prob);
            else return log(error_prob);
        }
    }
    else { // female or autosome
        if(obs_gen==true_gen) return log(1.0-error_prob);
        else return log(error_prob);
    }
}


double BC::step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female,
                IntegerVector cross_info)
{
    check_geno(gen_left, false, is_X_chr, is_female, cross_info);
    check_geno(gen_right, false, is_X_chr, is_female, cross_info);

    if(gen_left==gen_right) return log(1.0-rec_frac);
    else return log(rec_frac);
}

IntegerVector BC::possible_gen(bool is_X_chr, bool is_female,
                               IntegerVector cross_info)
{
    if(!is_X_chr || (is_X_chr && is_female)) {
        int vals[] = {AA,AB};
        IntegerVector result(vals, vals+2);
        return result;
    }
    else {
        int vals[] = {AY,BY};
        IntegerVector result(vals, vals+2);
        return result;
    }
}

int BC::ngen(bool is_X_chr)
{
    if(is_X_chr) return 4;
    return 2;
}


double BC::nrec(int gen_left, int gen_right,
                bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    check_geno(gen_left, false, is_X_chr, is_female, cross_info);
    check_geno(gen_right, false, is_X_chr, is_female, cross_info);

    if(gen_left == gen_right) return 0.0;
    else return 1.0;
}

double BC::est_rec_frac(NumericMatrix gamma, bool is_X_chr)
{

    int n_gen = gamma.rows();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0;
    for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

    double diagsum = 0.0;
    for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

    return(1.0 - diagsum/denom);
}

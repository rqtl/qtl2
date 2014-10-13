// backcross HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_bc.h"

enum gen {AA=1, AB=2, BB=3, AY=3, BY=4};

const bool BC::check_geno(const int gen, const bool is_observed_value,
                          const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info)
{
    if(is_observed_value && gen==0) return true;

    if(!is_x_chr || (is_x_chr && is_female)) {
        if(gen != AA && gen != AB) return false;
    }
    else { // male X chr
        if(is_observed_value) {
            if(gen != AA && gen != BB) return false;
        }
        else {
            if(gen != AY && gen != BY) return false;
        }
    }

    return true;
}

const double BC::init(const int true_gen,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return log(0.5);
}

const double BC::emit(const int obs_gen, const int true_gen, const double error_prob,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
        return 0.0; // missing or invalid

    if(is_x_chr && !is_female) { // X chr males different
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


const double BC::step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool& is_x_chr, const bool& is_female,
                      const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left==gen_right) return log(1.0-rec_frac);
    else return log(rec_frac);
}

const IntegerVector BC::possible_gen(const bool& is_x_chr, const bool& is_female,
                                     const IntegerVector& cross_info)
{
    if(!is_x_chr || (is_x_chr && is_female)) {
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

const int BC::ngen(const bool& is_x_chr)
{
    if(is_x_chr) return 4;
    return 2;
}


const double BC::nrec(const int gen_left, const int gen_right,
                      const bool& is_x_chr, const bool& is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left == gen_right) return 0.0;
    else return 1.0;
}

const double BC::est_rec_frac(const NumericMatrix& gamma, const bool& is_x_chr)
{

    int n_gen = gamma.rows();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0;
    for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

    double diagsum = 0.0;
    for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

    return 1.0 - diagsum/denom;
}

// intercross QTLCross class (for HMM)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_f2.h"

enum gen {AA=1, AB=2, BB=3, notA=5, notB=4,
          AAX=1, ABX=2, BAX=3, BBX=4, AY=5, BY=6};

const bool F2::check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==AA || gen==AB || gen==BB ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    if(is_x_chr) {
        bool forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(forward_direction && (gen==AAX || gen==ABX)) return true;
            if(!forward_direction && (gen==BAX || gen==BBX)) return true;
        }
        else if(gen==AY || gen==BY) return true;
    }
    else if(gen==AA || gen==AB || gen==BB) return true;

    return false; // otherwise a problem
}

const double F2::init(const int true_gen,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr) return log(0.5);
    else {
        if(true_gen==AB) return log(0.5);
        else return log(0.25);
    }
}

const double F2::emit(const int obs_gen, const int true_gen, const double error_prob,
                      const IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
       return 0.0; // missing or invalid

    if(is_x_chr) {
        if(is_female) {
            bool is_forward_direction = (cross_info[0] == 0);
            if(is_forward_direction) {
                if(true_gen==AAX) {
                    if(obs_gen==AA) return log(1.0-error_prob);
                    if(obs_gen==AB || obs_gen==notA) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
                if(true_gen==ABX) {
                    if(obs_gen==AB || obs_gen==notA) return log(1.0-error_prob);
                    if(obs_gen==AA) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
            }
            else {
                if(true_gen==BAX) {
                    if(obs_gen==AB || obs_gen==notB) return log(1.0-error_prob);
                    if(obs_gen==BB) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
                if(true_gen==BBX) {
                    if(obs_gen==BB) return log(1.0-error_prob);
                    if(obs_gen==AB || obs_gen==notB) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
            }
        }
        else { // males
            if(true_gen==AY) {
                if(obs_gen==AA || obs_gen==notB) return log(1.0-error_prob);
                if(obs_gen==BB || obs_gen==notA) return log(error_prob);
                return 0.0; // treat everything else as missing
            }
            if(true_gen==BY) {
                if(obs_gen==BB || obs_gen==notA) return log(1.0-error_prob);
                if(obs_gen==AA || obs_gen==notB) return log(error_prob);
                return 0.0; // treat everything else as missing
            }
        }
    }
    else { // autosome
        if(true_gen==AA) {
            if(obs_gen==AA) return log(1.0-error_prob);
            if(obs_gen==AB || obs_gen==BB) return log(error_prob/2.0);
            if(obs_gen==notB) return log(1.0-error_prob/2.0);
            if(obs_gen==notA) return log(error_prob);
        }
        if(true_gen==AB) {
            if(obs_gen==AB) return log(1.0-error_prob);
            if(obs_gen==AA || obs_gen==BB) return log(error_prob/2.0);
            if(obs_gen==notB || obs_gen==notA) return log(1.0-error_prob/2.0);
        }
        if(true_gen==BB) {
            if(obs_gen==BB) return log(1.0-error_prob);
            if(obs_gen==AB || obs_gen==AA) return log(error_prob/2.0);
            if(obs_gen==notA) return log(1.0-error_prob/2.0);
            if(obs_gen==notB) return log(error_prob);
        }
    }

    return NA_REAL; // shouldn't get here
}


const double F2::step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr) {
        if(gen_left == gen_right) return log(1.0-rec_frac);
        else return log(rec_frac);
    }
    else { // autosome
        switch(gen_left) {
        case AA:
            switch(gen_right) {
            case AA: return 2.0*log(1.0-rec_frac);
            case AB: return log(2.0)+log(1.0-rec_frac)+log(rec_frac);
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
            case AB: return log(2.0)+log(1.0-rec_frac)+log(rec_frac);
            case BB: return 2.0*log(1.0-rec_frac);
            }
        }
    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector F2::possible_gen(const bool is_x_chr, const bool is_female,
                                     const IntegerVector& cross_info)
{
    if(is_x_chr) {
        bool is_forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(is_forward_direction) {
                int vals[] = {AAX,ABX};
                IntegerVector result(vals, vals+2);
                return result;
            }
            else {
                int vals[] = {BAX,BBX};
                IntegerVector result(vals, vals+2);
                return result;
            }
        }
        else { // male
            int vals[] = {AY,BY};
            IntegerVector result(vals, vals+2);
            return result;
        }
    }
    else { // autosome
        int vals[] = {AA,AB,BB};
        IntegerVector result(vals, vals+3);
        return result;
    }
}

const int F2::ngen(const bool is_x_chr)
{
    if(is_x_chr) return 6;
    return 3;
}

const NumericMatrix F2::geno2allele_matrix(const bool is_x_chr)
{
    if(is_x_chr) // no conversion needed
        return NumericMatrix(0,0);

    NumericMatrix result(3,2);
    result(0,0) = 1.0;
    result(1,0) = result(1,1) = 0.5;
    result(2,1) = 1.0;

    return result;
}

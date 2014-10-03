// intercross HMM functions

#include "cross.h"

enum gen {GENO_NA=0, AA=1, AB=2, BB=3, notA=5, notB=4, AY=1, BY=3};

bool F2::check_geno(int gen, bool is_observed_value,
                    bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==GENO_NA || gen==AA || gen==AB || gen==BB || 
           gen==notA || gen==notB) return true;
        throw std::range_error("Invalid observed genotype");
    }

    if(is_X_chr) {
        bool forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(forward_direction && (gen==AA || gen==AB)) return true;
            if(!forward_direction && (gen==AB || gen==BB)) return true;
        }
        else if(gen==AY || gen==BY) return true;
    }
    else if(gen==AA || gen==AB || gen==BB) return true;

    throw std::range_error("Invalid genotype");
    return false; // can't get here
}

double F2::init(int true_gen,
                bool is_X_chr, bool is_female,
                IntegerVector cross_info)
{
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    if(is_X_chr) return log(0.5);
    else {
        if(true_gen==AB) return log(0.5);
        else return log(0.25);
    }
}

double F2::emit(int obs_gen, int true_gen, double error_prob,
                bool is_X_chr, bool is_female,
                IntegerVector cross_info)
{
    check_geno(obs_gen, true, is_X_chr, is_female, cross_info);
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    if(obs_gen==GENO_NA) return 0.0; // log(1.0)

    if(is_X_chr) {
        if(is_female) {
            bool is_forward_direction = (cross_info[0] == 0);
            if(is_forward_direction) {
                if(true_gen==AA) {
                    if(obs_gen==AA) return log(1.0-error_prob);
                    if(obs_gen==AB || obs_gen==notA) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
                if(true_gen==AB) {
                    if(obs_gen==AB || obs_gen==notA) return log(1.0-error_prob);
                    if(obs_gen==AA) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
            }
            else {
                if(true_gen==AB) {
                    if(obs_gen==AB || obs_gen==notB) return log(1.0-error_prob);
                    if(obs_gen==BB) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
                if(true_gen==BB) {
                    if(obs_gen==BB) return log(1.0-error_prob);
                    if(obs_gen==AB || obs_gen==notB) return log(error_prob);
                    return 0.0; // treat everything else as missing
                }
            }
        }
        else { // males
            if(true_gen==AY) {
                if(obs_gen==AY || obs_gen==notB) return log(1.0-error_prob);
                if(obs_gen==BY || obs_gen==notA) return log(error_prob);
                return 0.0; // treat everything else as missing
            }
            if(true_gen==BY) {
                if(obs_gen==BY || obs_gen==notA) return log(1.0-error_prob);
                if(obs_gen==AY || obs_gen==notB) return log(error_prob);
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

    throw std::range_error("invalid obs_gen or true_gen");
    return NA_REAL; // can't get here
}


double F2::step(int gen_left, int gen_right, double rec_frac,
                bool is_X_chr, bool is_female,
                IntegerVector cross_info)
{
    check_geno(gen_left, false, is_X_chr, is_female, cross_info);
    check_geno(gen_right, false, is_X_chr, is_female, cross_info);

    if(is_X_chr) {
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

    throw std::range_error("invalid genotypes.");
    return NA_REAL; // can't get here
}

IntegerVector F2::geno_index(bool is_X_chr, bool is_female,
                       IntegerVector cross_info)
{
    if(is_X_chr) {
        bool is_forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(is_forward_direction) {
                int vals[] = {AA-1,AB-1};
                IntegerVector result(vals, vals+2);
                return result;
            }
            else {
                int vals[] = {AB-1,BB-1};
                IntegerVector result(vals, vals+2);
                return result;
            }
        }
        else { // male
            int vals[] = {AA-1,BB-1};
            IntegerVector result(vals, vals+2);
            return result;
        }
    }
    else { // autosome
        int vals[] = {AA-1,AB-1,BB-1};
        IntegerVector result(vals, vals+3);
        return result;
    }
}

int F2::n_geno(bool is_X_chr)
{
    return 3;
}

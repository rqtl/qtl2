// intercross QTLCross class (for HMM)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_f2.h"
#include "r_message.h"

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
    #ifndef NDEBUG
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
    #ifndef NDEBUG
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
    #ifndef NDEBUG
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
                IntegerVector result = IntegerVector::create(AAX,ABX);
                return result;
            }
            else {
                IntegerVector result = IntegerVector::create(BAX,BBX);
                return result;
            }
        }
        else { // male
            IntegerVector result = IntegerVector::create(AY,BY);
            return result;
        }
    }
    else { // autosome
        IntegerVector result = IntegerVector::create(AA,AB,BB);
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
    if(is_x_chr) { // X chr
        NumericMatrix result(6,4);
        // female X
        result(0,0) = 1.0;
        result(1,0) = result(1,1) = 0.5;
        result(2,0) = result(2,1) = 0.5;
        result(3,1) = 1.0;

        // male X
        result(4,0) = result(5,1) = 1.0;

        return result;
    }
    else {
        NumericMatrix result(3,2);

        result(0,0) = 1.0;
        result(1,0) = result(1,1) = 0.5;
        result(2,1) = 1.0;

        return result;
    }
}

// check that sex conforms to expectation
const bool F2::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n = is_female.size();
    if(!any_x_chr) { // all autosomes
        if(n > 0) {
            result = true; // don't call this an error
            r_message("is_female included but not needed without X chromosome");
        }
    }
    else { // X chr included
        if(n == 0) {
            result = false;
            r_message("is_female not provided, but needed to handle X chromosome");
        }
        else {
            unsigned int n_missing = 0;
            for(unsigned int i=0; i<n; i++)
                if(is_female[i] == NA_LOGICAL) ++n_missing;
            if(n_missing > 0) {
                result = false;
                r_message("is_female contains missing values (it shouldn't)");
            }
        }
    }
    return result;
}

// check that cross_info conforms to expectation
const bool F2::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n_row = cross_info.rows();
    const unsigned int n_col = cross_info.cols();

    if(!any_x_chr) { // all autosomes
        if(n_col > 0) {
            result = true; // don't call this an error
            r_message("cross_info included but not needed without X chromosome");
        }
    }
    else { // X chr included
        if(n_col == 0) {
            result = false;
            r_message("cross_info not provided, but needed to handle X chromosome");
        }
        else if(n_col > 1) {
            result = false;
            r_message("cross_info has >1 columns, but should have just 1");
        }
        else {
            unsigned int n_missing = 0;
            for(unsigned int i=0; i<n_row; i++)
                if(cross_info[i] == NA_INTEGER) ++n_missing;
            if(n_missing > 0) {
                result = false;
                r_message("cross_info contains missing values (it shouldn't)");
            }

            unsigned int n_invalid = 0;
            for(unsigned int i=0; i<n_row; i++)
                if(cross_info[i] != NA_INTEGER && cross_info[i] != 0 && cross_info[i] != 1) ++n_invalid;
            if(n_invalid > 0) {
                result = false;
                r_message("cross_info contains invalid values; should be 0 or 1.");
            }
        }
    }
    return result;
}

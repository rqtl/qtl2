// backcross QTLCross class (for HMM)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_bc.h"
#include "r_message.h"

enum gen {AA=1, AB=2, BB=3, AY=3, BY=4};

const bool BC::check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
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
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return log(0.5);
}

const double BC::emit(const int obs_gen, const int true_gen, const double error_prob,
                      const IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const IntegerVector& cross_info)
{
    #ifndef NDEBUG
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
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left==gen_right) return log(1.0-rec_frac);
    else return log(rec_frac);
}

const IntegerVector BC::possible_gen(const bool is_x_chr, const bool is_female,
                                     const IntegerVector& cross_info)
{
    if(!is_x_chr || (is_x_chr && is_female)) {
        IntegerVector result = IntegerVector::create(AA,AB);
        return result;
    }
    else {
        IntegerVector result = IntegerVector::create(AY,BY);
        return result;
    }
}

const int BC::ngen(const bool is_x_chr)
{
    if(is_x_chr) return 4;
    return 2;
}


const double BC::nrec(const int gen_left, const int gen_right,
                      const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left == gen_right) return 0.0;
    else return 1.0;
}

// check that sex conforms to expectation
const bool BC::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n = is_female.size();
    if(!any_x_chr) { // all autosomes
        if(n > 0) {
            // not needed here, but don't call it an error
            result = true;
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

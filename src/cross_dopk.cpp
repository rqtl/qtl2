// phase-known Diversity Outcross QTLCross class (for HMM, in particular est.map)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_do.h"
#include "cross_dopk.h"
#include "cross_do_util.h"
#include "r_message.h"

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool DOPK::is_het(const int true_gen)
{
    IntegerVector alleles = decode_geno(true_gen);
    if(alleles[0] == alleles[1]) return false;
    return true;
}

// alleles -> integer 1, 2, ..., 36 (phase unknown case)
const int DOPK::encode_alleles(const int allele1, const int allele2)
{
    const int m = std::max(allele1, allele2);
    const int d = abs(allele1 - allele2);

    return (int)round(R::choose((double)(m+1), 2.0) - d);
}

// integer 1, 2, ..., 36 -> alleles (phase unknown case)
const IntegerVector DOPK::decode_geno(const int true_gen)
{
    const int n_geno = 36;
    #ifndef NDEBUG
    if(true_gen < 0 || true_gen > n_geno)
        throw std::range_error("genotype value not allowed");
    #endif

    IntegerVector result(2);

    int last_max = 0;
    for(int i=1; i<=n_geno; i++) {
        if(true_gen <= last_max+i) {
            result[1] = i;
            result[0] = true_gen - last_max;
            return(result);
        }
        last_max += i;
    }

    result[0] = NA_INTEGER;
    result[1] = NA_INTEGER;
    return result;
}



const bool DOPK::check_geno(const int gen, const bool is_observed_value,
                            const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // need to fill in this function
    return false;
}


const double DOPK::init(const int true_gen,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    // need to fill in this function
    return NA_REAL;
}

const double DOPK::emit(const int obs_gen, const int true_gen, const double error_prob,
                        const IntegerVector& founder_geno, const bool is_x_chr,
                        const bool is_female, const IntegerVector& cross_info)
{
    // need to fill in this function
    return NA_REAL;
}


const double DOPK::step(const int gen_left, const int gen_right, const double rec_frac,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    // need to fill in this function
    return NA_REAL;
}

const IntegerVector DOPK::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    int n_alleles = 8;
    int n_geno = 64;

    if(is_x_chr && !is_female) { // male X chromosome
        IntegerVector result(n_alleles);
        for(int i=0; i<n_alleles; i++)
            result[i] = n_geno+i+1;
        return result;
    }
    else { // autosome or female X
        IntegerVector result(n_geno);
        for(int i=0; i<n_geno; i++)
            result[i] = i+1;
        return result;
    }
}

const int DOPK::ngen(const bool is_x_chr)
{
    int n_alleles = 8;
    int n_geno = 64;

    if(is_x_chr) return n_geno+n_alleles;
    return n_geno;
}

const double DOPK::nrec(const int gen_left, const int gen_right,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    // need to fill in this function
    return NA_REAL;
}

const double DOPK::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                const IntegerMatrix& cross_info, const int n_gen)
{
    Rcpp::stop("est_map not yet available for Diversity Outcross");

    return NA_REAL;
}

const NumericMatrix DOPK::geno2allele_matrix(const bool is_x_chr)
{
    const int n_alleles = 8;
    const int n_geno = 64;

    if(is_x_chr) {
        NumericMatrix result(n_geno+n_alleles, n_alleles*2);
        // female X
        for(int trueg=0; trueg<n_geno; trueg++) {
            IntegerVector alleles = decode_geno(trueg+1);
            result(trueg,alleles[0]) += 0.5;
            result(trueg,alleles[1]) += 0.5;
        }
        // male X
        for(int trueg=0; trueg<n_geno; trueg++)
            result(trueg+n_geno, trueg+n_alleles) = 1.0;

        return result;
    }
    else { // autosome
        NumericMatrix result(n_geno,n_alleles);

        for(int trueg=0; trueg<n_geno; trueg++) {
            IntegerVector alleles = decode_geno(trueg+1);
            result(trueg,alleles[0]) += 0.5;
            result(trueg,alleles[1]) += 0.5;
        }

        return result;
    }
}

// check that sex conforms to expectation
const bool DOPK::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
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
const bool DOPK::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n_row = cross_info.rows();
    const unsigned int n_col = cross_info.cols();
    // one column with number of generations (needed no matter what; values should be >= 1)

    if(n_col == 0) {
        result = false;
        r_message("cross_info not provided, but should at least one column, with no. generations");
        return result;
    }

    unsigned int n_missing=0;
    unsigned int n_invalid=0;
    for(unsigned int i=0; i<n_row; i++) {
        if(cross_info[i] == NA_INTEGER) ++n_missing;
        else if(cross_info[i] < 1) ++n_invalid;
    }
    if(n_missing > 0) {
        result = false;
        r_message("cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("cross_info has invalid values; no. generations should be >= 1");
    }

    return result;
}

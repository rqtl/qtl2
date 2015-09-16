// Diversity Outcross QTLCross class (for HMM)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_do.h"
#include "cross_do_util.h"
#include "r_message.h"

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool DO::is_het(const int true_gen)
{
    IntegerVector alleles = DO::decode_geno(true_gen);
    if(alleles[0] == alleles[1]) return false;
    return true;
}

// alleles -> integer 1, 2, ..., 36 (phase unknown case)
const int DO::encode_alleles(const int allele1, const int allele2)
{
    int m = std::max(allele1, allele2);
    int d = abs(allele1 - allele2);

    return (int)round(R::choose((double)(m+1), 2.0) - d);
}

// integer 1, 2, ..., 36 -> alleles (phase unknown case)
const IntegerVector DO::decode_geno(const int true_gen)
{
    int n_alleles = 8;
    int n_geno = (int)round(R::choose((double)(n_alleles+1), 2.0));
    #ifdef DEBUG
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


const bool DO::check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    int n_alleles = 8;
    int n_geno = (int)round(R::choose((double)(n_alleles+1), 2.0));

    if(!is_x_chr || is_female) { // autosome or female X
        if(gen>= 1 && gen <= n_geno) return true;
    }
    else { // male X
        if(gen>=n_geno+1 && gen <=n_geno+n_alleles) return true;
    }

    return false; // otherwise a problem
}

const double DO::init(const int true_gen,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(!is_x_chr || is_female) { // autosome or female X
        if(is_het(true_gen)) return -log(32.0);
        else return -log(64.0);
    }
    else { // male X
        return -log(8.0);
    }
}

const double DO::emit(const int obs_gen, const int true_gen, const double error_prob,
                      const IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return NA_REAL; // shouldn't get here
}


const double DO::step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return NA_REAL; // shouldn't get here
}

const IntegerVector DO::possible_gen(const bool is_x_chr, const bool is_female,
                                     const IntegerVector& cross_info)
{
    int n_alleles = 8;
    int n_geno = (int)round(R::choose((double)(n_alleles+1), 2.0));

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

const int DO::ngen(const bool is_x_chr)
{
    int n_alleles = 8;
    int n_geno = (int)round(R::choose((double)(n_alleles+1), 2.0));

    if(is_x_chr) return n_geno+n_alleles;
    return n_geno;
}

const NumericMatrix DO::geno2allele_matrix(const bool is_x_chr)
{
    if(is_x_chr) // no conversion needed
        return NumericMatrix(0,0);

    NumericMatrix result(3,2);
    result(0,0) = 1.0;
    result(1,0) = result(1,1) = 0.5;
    result(2,1) = 1.0;

    return result;
}

// check that sex conforms to expectation
const bool DO::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
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
const bool DO::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n_row = cross_info.rows();
    const unsigned int n_col = cross_info.cols();
    // first column is number of generations (needed no matter what; values should be >= 2)
    // second column is 0=AxB, 1=BxA, 2=balanced (needed for X chromosome)

    if(n_col == 0) {
        result = false;
        r_message("cross_info not provided, but should at least have one column, with no. generations");
        return result;
    }

    unsigned int n_missing=0;
    unsigned int n_invalid=0;
    for(unsigned int i=0; i<n_row; i++) {
        if(cross_info[i] == NA_INTEGER) ++n_missing;
        else if(cross_info[i] < 2) ++n_invalid;
    }
    if(n_missing > 0) {
        result = false;
        r_message("1st column in cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("1st column in cross_info has invalid values; no. generations should be >= 2");
    }

    if(n_col == 1 && any_x_chr) {
        result = false;
        r_message("cross_info should have at two columns (no. generations and cross direction)");
    }

    if(n_col > 1) {
        if(n_col > 2) {
            result = false;
            r_message("cross_info should have no more than 2 columns (no. generations and cross direction)");
        }

        unsigned int n_missing = 0;
        unsigned int n_invalid = 0;
        for(unsigned int i=0; i<n_row; i++) {
            if(cross_info[i+n_row] == NA_INTEGER) ++n_missing;
            else if(cross_info[i+n_row] != 0 &&
                    cross_info[i+n_row] != 1 &&
                    cross_info[i+n_row] != 2) {
                ++n_invalid;
            }
        }
        if(n_missing > 0) {
            result = false;
            r_message("2nd column in cross_info contains missing values (it shouldn't)");
        }

        if(n_invalid > 0) {
            result = false;
            r_message("2nd column in cross_info contains invalid values; should be 0, 1, or 2.");
        }
    }
    return result;
}

// Heterogeneous Stock QTLCross class (for HMM)

#include "cross_hs.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_do_util.h"
#include "r_message.h"

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool HS::check_geno(const int gen, const bool is_observed_value,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_alleles = 8;
    const int n_geno = 36;

    if(!is_x_chr || is_female) { // autosome or female X
        if(gen>= 1 && gen <= n_geno) return true;
    }
    else { // male X
        if(gen>=n_geno+1 && gen <=n_geno+n_alleles) return true;
    }

    return false; // otherwise a problem
}

const double HS::init(const int true_gen,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(!is_x_chr || is_female) { // autosome or female X
        if(mpp_is_het(true_gen, 8, false)) return -log(32.0);
        else return -log(64.0);
    }
    else { // male X
        return -log(8.0);
    }
}

const double HS::emit(const int obs_gen, const int true_gen, const double error_prob,
                      const IntegerVector& founder_geno, const bool is_x_chr,
                      const bool is_female, const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif
    const int n_geno = 36;

    if(obs_gen==0) return 0.0; // missing

    if(!is_x_chr || is_female) { // autosome or female X
        const IntegerVector true_alleles = mpp_decode_geno(true_gen, 8, false);
        int f1 = founder_geno[true_alleles[0]-1];
        int f2 = founder_geno[true_alleles[1]-1];

        // treat founder hets as missing
        if(f1==2) f1 = 0;
        if(f2==2) f2 = 0;

        // neither founder alleles observed
        if(f1 == 0 && f2 == 0) return 0.0;

        // one founder allele observed
        if(f1 == 0 || f2 == 0) {

            switch(std::max(f1, f2)) {
            case H: return 0.0; // het compatible with either founder allele
            case A:
                switch(obs_gen) {
                case A: case notB: return log(1.0-error_prob);
                case B: case notA: return log(error_prob);
                case H: return 0.0;
                }
            case B:
                switch(obs_gen) {
                case B: case notA: return log(1.0-error_prob);
                case A: case notB: return log(error_prob);
                case H: return 0.0;
                }
            }
            return 0.0;
        }

        switch((f1+f2)/2) { // values 1, 2, 3
        case A:
            switch(obs_gen) {
            case A: return log(1.0-error_prob);
            case H: return log(error_prob/2.0);
            case B: return log(error_prob/2.0);
            case notA: return log(error_prob);
            case notB: return log(1.0-error_prob/2.0);
            }
        case H:
            switch(obs_gen) {
            case A: return log(error_prob/2.0);
            case H: return log(1.0-error_prob);
            case B: return log(error_prob/2.0);
            case notA: return log(1.0-error_prob/2.0);
            case notB: return log(1.0-error_prob/2.0);
            }
        case B:
            switch(obs_gen) {
            case B: return log(1.0-error_prob);
            case H: return log(error_prob/2.0);
            case A: return log(error_prob/2.0);
            case notB: return log(error_prob);
            case notA: return log(1.0-error_prob/2.0);
            }
        }
        return 0.0;
    }
    else { // male X
        const int founder_allele = founder_geno[(true_gen - n_geno) - 1];

        switch(founder_allele) {
        case A:
            switch(obs_gen) {
            case A: case notB: return log(1.0-error_prob);
            case B: case notA: return log(error_prob);
            }
        case B:
            switch(obs_gen) {
            case B: case notA: return log(1.0-error_prob);
            case A: case notB: return log(error_prob);
            }
        }
        return(0.0);
    }

    return NA_REAL; // shouldn't get here
}


const double HS::step(const int gen_left, const int gen_right, const double rec_frac,
                      const bool is_x_chr, const bool is_female,
                      const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    // can treat HS like DO where "preCC" founders all at generation 1
    const static IntegerVector precc_gen = IntegerVector::create(1);
    const static NumericVector precc_alpha =
        NumericVector::create(1.0);

    // no. generations for this mouse
    int n_gen = cross_info[0];

    if(is_x_chr) {
        if(is_female) { // female X
            return DOstep_femX(gen_left, gen_right, rec_frac, n_gen,
                               precc_gen, precc_alpha);
        }
        else { // male X
            return DOstep_malX(gen_left, gen_right, rec_frac, n_gen,
                               precc_gen, precc_alpha);
        }
    }
    else { // autosome
        return DOstep_auto(gen_left, gen_right, rec_frac, n_gen,
                           precc_gen, precc_alpha);
    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector HS::possible_gen(const bool is_x_chr, const bool is_female,
                                     const IntegerVector& cross_info)
{
    int n_alleles = 8;
    int n_geno = 36;

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

const int HS::ngen(const bool is_x_chr)
{
    int n_alleles = 8;
    int n_geno = 36;

    if(is_x_chr) return n_geno+n_alleles;
    return n_geno;
}

const int HS::nalleles()
{
    return 8;
}

const NumericMatrix HS::geno2allele_matrix(const bool is_x_chr)
{
    const int n_alleles = 8;
    const int n_geno = 36;

    if(is_x_chr) {
        NumericMatrix result(n_geno+n_alleles, n_alleles);
        // female X
        for(int trueg=0; trueg<n_geno; trueg++) {
            IntegerVector alleles = mpp_decode_geno(trueg+1, 8, false);
            result(trueg,alleles[0]-1) += 0.5;
            result(trueg,alleles[1]-1) += 0.5;
        }
        // male X
        for(int trueg=0; trueg<n_alleles; trueg++)
            result(trueg+n_geno, trueg) = 1.0;

        return result;
    }
    else { // autosome
        NumericMatrix result(n_geno,n_alleles);

        for(int trueg=0; trueg<n_geno; trueg++) {
            IntegerVector alleles = mpp_decode_geno(trueg+1, 8, false);
            result(trueg,alleles[0]-1) += 0.5;
            result(trueg,alleles[1]-1) += 0.5;
        }

        return result;
    }
}

// check that sex conforms to expectation
const bool HS::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
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

// check that cross_info conforms to expectation
const bool HS::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
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

// check that founder genotype data has correct no. founders and markers
const bool HS::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != 8) {
        result = false;
        r_message("founder_geno should have 8 founders");
    }

    return result;
}

// check that founder genotype data has correct values
const bool HS::check_founder_geno_values(const IntegerMatrix& founder_geno)
{
    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    for(int f=0; f<fg_f; f++) {
        for(int mar=0; mar<fg_mar; mar++) {
            int fg = founder_geno(f,mar);
            if(fg != 0 && fg != 1 && fg != 3) {
                // at least one invalid value
                r_message("founder_geno contains invalid values; should be in {0, 1, 3}");
                return false;
            }
        }
    }

    return true;
}

const bool HS::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> HS::geno_names(const std::vector<std::string> alleles,
                                              const bool is_x_chr)
{
    if(alleles.size() < 8)
        throw std::range_error("alleles must have length 8");

    if(is_x_chr) {
        std::vector<std::string> result(44);
        for(int i=0; i<36; i++) {
            IntegerVector allele_int = mpp_decode_geno(i+1, 8, false);
            result[i] = alleles[allele_int[0]-1] + alleles[allele_int[1]-1];
        }
        for(int i=0; i<8; i++) {
            result[36+i] = alleles[i] + "Y";
        }

        return result;
    }
    else {
        std::vector<std::string> result(36);
        for(int i=0; i<36; i++) {
            IntegerVector allele_int = mpp_decode_geno(i+1, 8, false);
            result[i] = alleles[allele_int[0]-1] + alleles[allele_int[1]-1];
        }
        return result;
    }
}

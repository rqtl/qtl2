// Diversity Outcross F1 (in cross with another inbred strain) QTLCross class (for HMM)

#include "cross_dof1.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_do_util.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen_dof1 {A=1, H=2, B=3, notA=5, notB=4};

const bool DOF1::check_geno(const int gen, const bool is_observed_value,
                            const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_alleles = 8;
    const int n_geno = 8;

    if(gen>= 1 && gen <= n_geno) return true;

    return false; // otherwise a problem
}

const double DOF1::init(const int true_gen,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return -log(8.0);
}

const double DOF1::emit(const int obs_gen, const int true_gen, const double error_prob,
                        const IntegerVector& founder_geno, const bool is_x_chr,
                        const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0) return 0.0; // missing

    int f1 = founder_geno[true_gen-1];

    int f2;
    if(is_x_chr && !is_female) {
        f2 = f1; // assuming female DO x male 9th strain, so male is hemizygous DO chr
    } else {
        // autosome or female X chr
        f2 = founder_geno[8]; // 9th founder is the other strain
    }

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


const double DOF1::step(const int gen_left, const int gen_right, const double rec_frac,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    // info about preCC progenitors
    const static IntegerVector precc_gen = IntegerVector::create(4,5,6,7,8,9,10,11,12);
    const static NumericVector precc_alpha =
        NumericVector::create(21.0/144.0, 64.0/144.0, 24.0/144.0, 10.0/144.0, 5.0/144.0,
                               9.0/144.0,  5.0/144.0,  3.0/144.0,  3.0/144.0);

    // no. generations for this mouse
    int n_gen = cross_info[0];

    double r;
    if(is_x_chr) {
        if(is_female) { // female X
            r = DOrec_femX(rec_frac, n_gen, precc_gen, precc_alpha);
        }
        else { // male X
            r = DOrec_malX(rec_frac, n_gen, precc_gen, precc_alpha);
        }
    }
    else { // autosome
        r = DOrec_auto(rec_frac, n_gen, precc_gen, precc_alpha);
    }

    if(gen_left == gen_right) return(log(1.0-r));
    else return(log(r)-log(7.0));
}

const IntegerVector DOF1::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    int n_geno = 8;

    IntegerVector result(n_geno);
    for(int i=0; i<n_geno; i++)
        result[i] = i+1;
    return result;

}

const int DOF1::ngen(const bool is_x_chr)
{
    return 8;
}

const int DOF1::nalleles()
{
    return 8;
}

const NumericMatrix DOF1::geno2allele_matrix(const bool is_x_chr)
{
    const int n_geno = 8;

    NumericMatrix result(n_geno,n_geno);

    for(int trueg=0; trueg<n_geno; trueg++)
        result(trueg,trueg) = 1.0;

    return result;
}

// check that sex conforms to expectation
const bool DOF1::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const int n = is_female.size();
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
            int n_missing = 0;
            for(int i=0; i<n; i++)
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
const bool DOF1::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();
    // one column with number of generations (needed no matter what; values should be >= 1)

    if(n_col == 0) {
        result = false;
        r_message("cross_info should have at least one column, with no. generations");
        return result;
    }

    int n_missing=0;
    int n_invalid=0;
    for(int i=0; i<n_row; i++) {
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
const bool DOF1::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != 9) {
        result = false;
        r_message("founder_geno should have 9 founders (the CC founders plus the 9th strain to which the DO individuals were crossed)");
    }

    return result;
}

// check that founder genotype data has correct values
const bool DOF1::check_founder_geno_values(const IntegerMatrix& founder_geno)
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

const bool DOF1::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> DOF1::geno_names(const std::vector<std::string> alleles,
                                                const bool is_x_chr)
{
    int n_alleles = alleles.size();

    std::vector<std::string> result(n_alleles);

    for(int i=0; i<n_alleles; i++)
        result[i] = alleles[i];

    return result;
}

const double DOF1::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                const IntegerMatrix& cross_info, const int n_gen)
{
    Rcpp::stop("est_map not yet available for Diversity Outcross F1");

    return NA_REAL;
}

// not implemented
const NumericVector DOF1::est_map2(const IntegerMatrix& genotypes,
                                   const IntegerMatrix& founder_geno,
                                   const bool is_X_chr,
                                   const LogicalVector& is_female,
                                   const IntegerMatrix& cross_info,
                                   const IntegerVector& cross_group,
                                   const IntegerVector& unique_cross_group,
                                   const NumericVector& rec_frac,
                                   const double error_prob,
                                   const int max_iterations,
                                   const double tol,
                                   const bool verbose)
{
    Rcpp::stop("est_map not yet implemented for DO F1s.");

    // return vector of NAs
    const int n_rf = rec_frac.size();
    NumericVector result(n_rf);
    for(int i=0; i<n_rf; i++) result[i] = NA_REAL;
    return result ;
}

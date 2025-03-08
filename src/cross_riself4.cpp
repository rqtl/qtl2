// 4-way RIL by selfing QTLCross class (for HMM)

#include "cross_riself4.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_do_util.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen_riself4 {A=1, H=2, B=3, notA=5, notB=4};

const bool RISELF4::check_geno(const int gen, const bool is_observed_value,
                                const bool is_x_chr, const bool is_female,
                                const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_geno = 4;

    if(gen>= 1 && gen <= n_geno) return true;

    return false; // otherwise a problem
}

const double RISELF4::init(const int true_gen,
                            const bool is_x_chr, const bool is_female,
                            const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return -log(4.0);
}

const double RISELF4::emit(const int obs_gen, const int true_gen, const double error_prob,
                            const IntegerVector& founder_geno, const bool is_x_chr,
                            const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0) return 0.0; // missing

    int f = founder_geno[true_gen-1]; // founder allele
    if(f!=1 && f!=3) return 0.0;      // founder missing -> no information

    if(f == obs_gen) return log(1.0 - error_prob);

    return log(error_prob); // genotyping error
}


const double RISELF4::step(const int gen_left, const int gen_right, const double rec_frac,
                            const bool is_x_chr, const bool is_female,
                            const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    // equations are from Broman (2005) Genetics 169:1133-1146
    //    doi:10.1534/genetics.104.035212
    //    see equation in right column on page 1135
    //    (need to multiply by 4 to get conditional probabilities)
    if(gen_left != gen_right)
        return log(rec_frac) - log(1.0 + 2.0*rec_frac);
    else
        return log(1.0 - rec_frac) - log(1.0 + 2.0*rec_frac);
}

const IntegerVector RISELF4::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    int n_geno = 4;
    IntegerVector result(n_geno);

    for(int i=0; i<n_geno; i++) result[i] = i+1;
    return result;
}

const int RISELF4::ngen(const bool is_x_chr)
{
    return 4;
}

const int RISELF4::nalleles()
{
    return 4;
}


// check that cross_info conforms to expectation
const bool RISELF4::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();
    // 4 columns with order of cross

    if(n_col != 4) {
        result = false;
        r_message("cross_info should have 4 columns, indicating the order of the cross");
        return result;
    }

    int n_missing=0;
    int n_invalid=0;
    for(int i=0; i<n_row; i++) {
        for(int j=0; j<n_col; j++) {
            if(cross_info(i,j) == NA_INTEGER) ++n_missing;
            else if(cross_info(i,j) < 1 || cross_info(i,j)>n_col) ++n_invalid;
        }
        // count values 1..ncol
        IntegerVector counts(n_col);
        for(int j=0; j<n_col; j++) counts[j] = 0; // zero counts
        for(int j=0; j<n_col; j++) {
            if(cross_info(i,j) >= 1 && cross_info(i,j)<=n_col) // ignore if out of range
                ++counts[cross_info(i,j)-1]; // count values
        }
        for(int j=0; j<n_col; j++) {
            if(counts[j] != 1) n_invalid += abs(counts[j] - 1);
        }
    }
    if(n_missing > 0) {
        result = false;
        r_message("cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("cross_info has invalid values; each row should be permutation of {1, 2, ..., 4}");
    }

    return result;
}


// check that founder genotype data has correct no. founders and markers
const bool RISELF4::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != 4) {
        result = false;
        r_message("founder_geno should have 4 founders");
    }

    return result;
}

// check that founder genotype data has correct values
const bool RISELF4::check_founder_geno_values(const IntegerMatrix& founder_geno)
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

const bool RISELF4::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> RISELF4::geno_names(const std::vector<std::string> alleles,
                                                const bool is_x_chr)
{
    if(alleles.size() < 4)
        throw std::range_error("alleles must have length 4");

    const int n_alleles = 4;

    std::vector<std::string> result(n_alleles);

    for(int i=0; i<n_alleles; i++)
        result[i] = alleles[i] + alleles[i];

    return result;
}


const int RISELF4::nrec(const int gen_left, const int gen_right,
                         const bool is_x_chr, const bool is_female,
                         const Rcpp::IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(gen_left == gen_right) return 0;
    else return 1;
}

const double RISELF4::est_rec_frac(const Rcpp::NumericVector& gamma, const bool is_x_chr,
                                    const Rcpp::IntegerMatrix& cross_info, const int n_gen)
{
    double R = QTLCross::est_rec_frac(gamma, is_x_chr, cross_info, n_gen);

    // inverse of R = 3r/(1+2r)
    return R/(3.0 - 2.0*R);
}

// check whether X chr can be handled
const bool RISELF4::check_handle_x_chr(const bool any_x_chr)
{
    if(any_x_chr) {
        r_message("X chr ignored for RIL by selfing.");
        return false;
    }

    return true; // most crosses can handle the X chr
}

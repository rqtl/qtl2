// QTLCross class for general k-way recombinant inbred lines or doubled haploids

#include "cross_genril.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool GENRIL::check_geno(const int gen, const bool is_observed_value,
                           const bool is_x_chr, const bool is_female,
                           const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    if(gen>= 1 && gen <= this->n_founders) return true;

    return false; // otherwise a problem
}

const double GENRIL::init(const int true_gen,
                       const bool is_x_chr, const bool is_female,
                       const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    // cross_info = num_generations followed by un-scaled frequences, as integers
    // so first sum up the values
    int denom = 0;
    for(int i=1; i<=this->n_founders; i++) denom += cross_info[i];

    // log frequency - log denominator
    return log((double)cross_info[true_gen]) - log((double)denom);
}

const double GENRIL::emit(const int obs_gen, const int true_gen, const double error_prob,
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


const double GENRIL::step(const int gen_left, const int gen_right, const double rec_frac,
                          const bool is_x_chr, const bool is_female,
                          const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    return step_genchr(gen_left, gen_right, rec_frac, is_x_chr, cross_info, this->n_founders);

}

const IntegerVector GENRIL::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    int n_geno = this->n_founders;
    IntegerVector result(n_geno);

    for(int i=0; i<n_geno; i++) result[i] = i+1;
    return result;
}

const int GENRIL::ngen(const bool is_x_chr)
{
    return this->n_founders;
}

const int GENRIL::nalleles()
{
    return this->n_founders;
}


// check that cross_info conforms to expectation
const bool GENRIL::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();

    // number of generations of outbreeding followed by rel frequency of each founder, as integers
    if(n_col != 1 + this->n_founders) {
        result = false;
        r_message("cross_info should have (1 + n_founders) columns: no. generations + rel freq of founders, as integers");
        return result;
    }

    int n_missing=0;
    int n_invalid=0;
    int n_zerosum=0;
    for(int i=0; i<n_row; i++) {
        if(cross_info(i,0) == NA_INTEGER) ++n_missing;
        else if(cross_info(i,0) < 2) ++n_invalid;
        int sum_alpha=0;
        for(int j=1; j<=this->n_founders; j++) {
            if(cross_info(i,j) == NA_INTEGER) ++n_missing;
            else if(cross_info(i,j) < 0) ++n_invalid;
            sum_alpha += cross_info(i,j);
        }
        if(sum_alpha == 0) n_zerosum++;
    }
    if(n_missing > 0) {
        result = false;
        r_message("cross_info has missing values (it shouldn't)");
    }
    if(n_invalid > 0) {
        result = false;
        r_message("cross_info has invalid values; no. gen should be >= 2 and rel freq should be >= 0");
    }
    if(n_zerosum > 0) {
        result = false;
        r_message("cross_info has invalid rows; rel freq should have positive sums");
    }

    return result;
}


// check that founder genotype data has correct no. founders and markers
const bool GENRIL::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
{
    bool result=true;

    const int fg_mar = founder_geno.cols();
    const int fg_f   = founder_geno.rows();

    if(fg_mar != n_markers) {
        result = false;
        r_message("founder_geno has incorrect number of markers");
    }

    if(fg_f != this->n_founders) {
        result = false;
        r_message("no. columns in founder_geno doesn't match no. founders");
    }

    return result;
}

// check that founder genotype data has correct values
const bool GENRIL::check_founder_geno_values(const IntegerMatrix& founder_geno)
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

const bool GENRIL::need_founder_geno()
{
    return true;
}

// geno_names from allele names
const std::vector<std::string> GENRIL::geno_names(const std::vector<std::string> alleles,
                                                const bool is_x_chr)
{
    if(alleles.size() < (unsigned)(this->n_founders))
        throw std::range_error("alleles must have length n_founders");

    const int n_alleles = this->n_founders;

    std::vector<std::string> result(n_alleles);

    for(int i=0; i<n_alleles; i++)
        result[i] = alleles[i] + alleles[i];

    return result;
}


const int GENRIL::nrec(const int gen_left, const int gen_right,
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


// check whether X chr can be handled
const bool GENRIL::check_handle_x_chr(const bool any_x_chr)
{
    return true; // most crosses can handle the X chr
}

// tailored est_map that pre-calculates transition matrices, etc
const NumericVector GENRIL::est_map2(const IntegerMatrix& genotypes,
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
    Rcpp::stop("est_map not yet implemented for general RIL.");

    // return vector of NAs
    const int n_rf = rec_frac.size();
    NumericVector result(n_rf);
    for(int i=0; i<n_rf; i++) result[i] = NA_REAL;
    return result ;
}

// step for general chromosome, used for both general RIL and general AIL
// TODO: at present, returns the same thing whether autosome or X chromosome
// TODO:    - could use X-chr-specific results
const double step_genchr(const int gen_left, const int gen_right, const double rec_frac,
                         const bool is_x_chr, const IntegerVector& cross_info, const int n_founders)
{
    #ifndef RQTL2_NODEBUG
    if(gen_left < 1 || gen_left > n_founders ||
       gen_right < 1 || gen_right > n_founders) {
        throw std::range_error("genotype value not allowed");
    }
    #endif

    const int k = cross_info[0]; // number of generations

    // sum of relative frequencies
    int denom=0.0;
    for(int i=0; i<n_founders; i++) denom += cross_info[i+1];

    if(gen_left == gen_right)
        return log((double)cross_info[gen_left] + pow(1.0-rec_frac, (double)k) * (denom - cross_info[gen_left])) -
            log((double)denom);
    else
        return log((double)cross_info[gen_right]) - log((double)denom) +
            log(1.0 - pow(1.0 - rec_frac, (double)k));

}

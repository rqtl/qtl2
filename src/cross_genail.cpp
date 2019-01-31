// QTLCross class for general k-way advanced intercross lines

#include "cross_genail.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"
#include "cross_genril.h" // defines step_genchr()
#include "hmm_util.h" // defines addlog()
#include "r_message.h" // defines RQTL2_NODEBUG and r_message()

enum gen {A=1, H=2, B=3, notA=5, notB=4};

const bool GENAIL::check_geno(const int gen, const bool is_observed_value,
                           const bool is_x_chr, const bool is_female,
                           const IntegerVector& cross_info)
{
    // allow any value 0-5 for observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    const int n_auto_geno = this->ngen(false);

    if(is_x_chr && !is_female) { // X chromosome, male
        if(gen>= n_auto_geno + 1 && gen <= n_auto_geno + this->n_founders) return true;
    }
    else {
        if(gen>= 1 && gen <= n_auto_geno) return true;
    }

    return false; // otherwise a problem
}

const double GENAIL::init(const int true_gen,
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

    const int n_auto_geno = this->ngen(false);

    if(is_x_chr && !is_female) { // male X chr
        // log frequency - log denominator
        return log((double)cross_info[true_gen-n_auto_geno]) - log((double)denom);
    }
    else { // autosome or female X
        const IntegerVector alleles = mpp_decode_geno(true_gen, this->n_founders, false);

        if(mpp_is_het(true_gen, this->n_founders, false)) {
            return log(2.0) + log((double)cross_info[alleles[0]]) +
                log((double)cross_info[alleles[1]]) - 2.0*log((double)denom);
        }
        else { // homozygote
            return (log((double)cross_info[alleles[0]]) - log((double)denom))*2.0;
        }
    }
}

// this basically follows the DO case
const double GENAIL::emit(const int obs_gen, const int true_gen, const double error_prob,
                       const IntegerVector& founder_geno, const bool is_x_chr,
                       const bool is_female, const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0) return 0.0; // missing

    const int n_auto_geno = this->ngen(false);

    if(!is_x_chr || is_female) { // autosome or female X
        const IntegerVector true_alleles = mpp_decode_geno(true_gen, this->n_founders, false);
        int f1 = founder_geno[true_alleles[0]-1];
        int f2 = founder_geno[true_alleles[1]-1];

        // treat founder hets as missing
        if(f1==2) f1 = 0;
        if(f2==2) f2 = 0;

        // neither founder alleles observed
        if(f1==0 && f2==0) return 0.0;

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
        const int founder_allele = founder_geno[(true_gen - n_auto_geno) - 1];

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
}


const double GENAIL::step(const int gen_left, const int gen_right, const double rec_frac,
                            const bool is_x_chr, const bool is_female,
                            const IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr && !is_female) { // male X chr

        const int n_auto_geno = this->ngen(false);

        if(gen_left < n_auto_geno || gen_right < n_auto_geno) {
            throw std::range_error("genotype value not allowed");
        }

        return step_genchr(gen_left - n_auto_geno, gen_right - n_auto_geno, rec_frac,
                           is_x_chr, cross_info, this->n_founders);
    }
    else { // autosome or female X

        const IntegerVector leftg = mpp_decode_geno(gen_left, this->n_founders, false);
        const IntegerVector rightg = mpp_decode_geno(gen_right, this->n_founders, false);

        const int left1 = leftg[0];
        const int left2 = leftg[1];
        const int right1 = rightg[0];
        const int right2 = rightg[1];

        if(left1 == left2) { // AA ->
            if(right1 == right2) { // AA -> AA (or AA -> BB)
                return 2.0*step_genchr(left1, right1, rec_frac, is_x_chr,
                                       cross_info, this->n_founders);
            }
            else { // AA -> AB (or AA -> BC)
                return step_genchr(left1, right1, rec_frac, is_x_chr,
                                   cross_info, this->n_founders) +
                    step_genchr(left1, right2, rec_frac, is_x_chr,
                                cross_info, this->n_founders);
            }
        }
        else { // AB ->
            if(right1 == right2) { // AB -> AA (or AB -> CC)
                return step_genchr(left1, right1, rec_frac, is_x_chr,
                                   cross_info, this->n_founders) +
                    step_genchr(left2, right1, rec_frac, is_x_chr,
                                cross_info, this->n_founders);
            }
            else { // AB -> CD (or AB -> AB or AB -> AC)
                // 1-1, 2-2 or 1-2, 2-1
                return addlog(step_genchr(left1, right1, rec_frac, is_x_chr,
                                          cross_info, this->n_founders),
                              step_genchr(left2, right2, rec_frac, is_x_chr,
                                          cross_info, this->n_founders)) +
                    addlog(step_genchr(left1, right2, rec_frac, is_x_chr,
                                       cross_info, this->n_founders),
                           step_genchr(left1, right2, rec_frac, is_x_chr,
                                       cross_info, this->n_founders));
            }
        }
    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector GENAIL::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    const int n_auto_geno = this->ngen(false);

    if(is_x_chr && !is_female) { // X chr male
        int n_geno = this->n_founders;
        IntegerVector result(n_geno);

        for(int i=0; i<n_geno; i++) result[i] = i+n_auto_geno+1;
        return result;
    }
    else { // autosome or X chr female
        IntegerVector result(n_auto_geno);

        for(int i=0; i<n_auto_geno; i++) result[i] = i+1;
        return result;
    }
}

const int GENAIL::ngen(const bool is_x_chr)
{
    const int n_auto_geno = this->n_founders + this->n_founders*(this->n_founders-1)/2;

    if(is_x_chr) {
        return n_auto_geno + this->n_founders;
    }
    else {
        return n_auto_geno;
    }
}

const int GENAIL::nalleles()
{
    return this->n_founders;
}

// check that cross_info conforms to expectation
const bool GENAIL::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const int n_row = cross_info.rows();
    const int n_col = cross_info.cols();
    // single column with the number of generations

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
const bool GENAIL::check_founder_geno_size(const IntegerMatrix& founder_geno, const int n_markers)
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
const bool GENAIL::check_founder_geno_values(const IntegerMatrix& founder_geno)
{
    const int n_f = this->n_founders;
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

const bool GENAIL::need_founder_geno()
{
    return true;
}


// geno_names from allele names
const std::vector<std::string> GENAIL::geno_names(const std::vector<std::string> alleles,
                                                const bool is_x_chr)
{
    if(alleles.size() != this->n_founders)
        throw std::range_error("alleles must have length n_founders");

    return mpp_geno_names(alleles, is_x_chr);
}


const int GENAIL::nrec(const int gen_left, const int gen_right,
                         const bool is_x_chr, const bool is_female,
                         const Rcpp::IntegerVector& cross_info)
{
    #ifndef RQTL2_NODEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const int n_auto_geno = this->ngen(false); // number of autosomal genotypes

    if(is_x_chr && gen_left > n_auto_geno && gen_right > n_auto_geno) { // male X chr
        if(gen_left == gen_right) return(0);
        else return(1);
    }

    Rcpp::IntegerVector a_left = mpp_decode_geno(gen_left, this->n_founders, false);
    Rcpp::IntegerVector a_right = mpp_decode_geno(gen_right, this->n_founders, false);

    if(a_left[0] == a_right[0]) {
        if(a_left[1] == a_right[1]) return(0);
        else return(1);
    }
    else if(a_left[0] == a_right[1]) {
        if(a_left[1] == a_right[0]) return(0);
        else return(1);
    }
    else if(a_left[1] == a_right[0]) {
        return(1);
    }
    else if(a_left[1] == a_right[1]) {
        return(1);
    }
    else return(2);
}


// check whether X chr can be handled
const bool GENAIL::check_handle_x_chr(const bool any_x_chr)
{
    return true; // most crosses can handle the X chr
}

// tailored est_map that pre-calculates transition matrices, etc
const NumericVector GENAIL::est_map2(const IntegerMatrix& genotypes,
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

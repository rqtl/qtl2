// phase-known intercross QTLCross class (for HMM, in particular est.map)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_f2pk.h"

enum gen {AA=1, AB=2, BA=3, BB=4,
          A=1, H=2, B=3, notB=4, notA=5,
          AY=5, BY=6};

const bool F2PK::check_geno(const int gen, const bool is_observed_value,
                            const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 or observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    if(is_x_chr) {
        bool forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(forward_direction && (gen==AA || gen==AB)) return true;
            if(!forward_direction && (gen==BA || gen==BB)) return true;
        }
        else if(gen==AY || gen==BY) return true;
    }
    else { // autosome
        if(gen==AA || gen==AB || gen==BA || gen==BB) return true;
    }

    return false; // invalid
}


const double F2PK::init(const int true_gen,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr) return log(0.5);
    else return log(0.25);
}

const double F2PK::emit(const int obs_gen, const int true_gen, const double error_prob,
                        const IntegerVector& founder_geno, const bool is_x_chr,
                        const bool is_female, const IntegerVector& cross_info)
{

    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
        return 0.0; // missing or invalid

    if(is_x_chr) {
        if(is_female) {
            bool is_forward_direction = (cross_info[0] == 0);
            if(is_forward_direction) {
                if(true_gen==AA) {
                    if(obs_gen==A) return log(1.0-error_prob);
                    if(obs_gen==H || obs_gen==notA) return log(error_prob);
                    return 0.0; // treat everything else as NA
                }
                if(true_gen==AB) {
                    if(obs_gen==H || obs_gen==notA) return log(1.0-error_prob);
                    if(obs_gen==A) return log(error_prob);
                    return 0.0; // treat everything else as NA
                }
            }
            else {
                if(true_gen==BA) {
                    if(obs_gen==H || obs_gen==notB) return log(1.0-error_prob);
                    if(obs_gen==B) return log(error_prob);
                    return 0.0; // treat everything else as NA
                }
                if(true_gen==BB) {
                    if(obs_gen==B) return log(1.0-error_prob);
                    if(obs_gen==H || obs_gen==notB) return log(error_prob);
                    return 0.0; // treat everything else as NA
                }
            }
        }
        else { // males
            if(true_gen==AY) {
                if(obs_gen==A || obs_gen==notB) return log(1.0-error_prob);
                if(obs_gen==B || obs_gen==notA) return log(error_prob);
                return 0.0; // treat everything else as NA
            }
            else { // BY
                if(obs_gen==B || obs_gen==notA) return log(1.0-error_prob);
                if(obs_gen==A || obs_gen==notB) return log(error_prob);
                return 0.0; // treat everything else as NA
            }
        }
    }
    else { // autosome
        if(true_gen==AA) {
            if(obs_gen==A) return log(1.0-error_prob);
            if(obs_gen==H || obs_gen==B) return log(error_prob/2.0);
            if(obs_gen==notB) return log(1.0-error_prob/2.0);
            if(obs_gen==notA) return log(error_prob);
        }
        else if(true_gen == AB || true_gen==BA) {
            if(obs_gen==H) return log(1.0-error_prob);
            if(obs_gen==A || obs_gen==B) return log(error_prob/2.0);
            if(obs_gen==notB || obs_gen==notA) return log(1.0-error_prob/2.0);
        }
        else { // BB
            if(obs_gen==B) return log(1.0-error_prob);
            if(obs_gen==H || obs_gen==A) return log(error_prob/2.0);
            if(obs_gen==notA) return log(1.0-error_prob/2.0);
            if(obs_gen==notB) return log(error_prob);
        }
    }

    return NA_REAL; // shouldn't get here
}


const double F2PK::step(const int gen_left, const int gen_right, const double rec_frac,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifdef DEBUG
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
            case AB: case BA: return log(1.0-rec_frac)+log(rec_frac);
            case BB: return 2.0*log(rec_frac);
            }
        case AB:
            switch(gen_right) {
            case AA: case BB: return log(rec_frac)+log(1.0-rec_frac);
            case AB: return log((1.0-rec_frac)*(1.0-rec_frac));
            case BA: return log(rec_frac*rec_frac);
            }
        case BA:
            switch(gen_right) {
            case AA: case BB: return log(rec_frac)+log(1.0-rec_frac);
            case AB: return log(rec_frac*rec_frac);
            case BA: return log((1.0-rec_frac)*(1.0-rec_frac));
            }
        case BB:
            switch(gen_right) {
            case AA: return 2.0*log(rec_frac);
            case AB: case BA: return log(1.0-rec_frac)+log(rec_frac);
            case BB: return 2.0*log(1.0-rec_frac);
            }
        }
    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector F2PK::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    if(is_x_chr) {
        bool is_forward_direction = (cross_info[0]==0);
        if(is_female) {
            if(is_forward_direction) {
                int vals[] = {AA,AB};
                IntegerVector result(vals, vals+2);
                return result;
            }
            else {
                int vals[] = {BA,BB};
                IntegerVector result(vals, vals+2);
                return result;
            }
        }
        else { // male
            int vals[] = {AY,BY};
            IntegerVector result(vals, vals+2);
            return result;
        }
    }
    else { // autosome
        int vals[] = {AA,AB,BA,BB};
        IntegerVector result(vals, vals+4);
        return result;
    }
}

const int F2PK::ngen(const bool is_x_chr)
{
    if(is_x_chr) return 6;
    return 4;
}

const double F2PK::nrec(const int gen_left, const int gen_right,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr) {
        if(gen_left == gen_right) return 0.0;
        else return 1.0;
    }
    else { // autosome
        switch(gen_left) {
        case AA:
            switch(gen_right) {
            case AA: return 0.0;
            case AB: case BA: return 0.5;
            case BB: return 1.0;
            }
        case AB:
            switch(gen_right) {
            case AA: case BB: return 0.5;
            case AB: return 0.0;
            case BA: return 1.0;
            }
        case BA:
            switch(gen_right) {
            case AA: case BB: return 0.5;
            case BA: return 0.0;
            case AB: return 1.0;
            }
        case BB:
            switch(gen_right) {
            case AA: return 1.0;
            case AB: case BA: return 0.5;
            case BB: return 0.0;
            }
        }
    }

    return NA_REAL; // shouldn't get here
}

const double F2PK::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                const IntegerMatrix& cross_info, const int n_gen)
{

    if(is_x_chr)
        return QTLCross::est_rec_frac(gamma, is_x_chr, cross_info, n_gen);

    int n_ind = cross_info.cols();
    int n_gen_sq = n_gen*n_gen;

    // get counts
    NumericMatrix num_rec(n_gen, n_gen);
    IntegerVector empty(0);

    for(int il=0; il<n_gen; il++) {
        for(int ir=0; ir<n_gen; ir++)
            num_rec(il,ir) = nrec(il+1, ir+1, false, false, empty);
    }

    double numerator=0.0;

    for(int ind=0, offset=0; ind<n_ind; ind++, offset += n_gen_sq)
        for(int il=0; il<n_gen; il++)
            for(int ir=0; ir<n_gen; ir++)
                numerator += gamma[offset + ir*n_gen + il]*num_rec(il,ir);

    return numerator/(double)n_ind;
}

const NumericMatrix F2PK::geno2allele_matrix(const bool is_x_chr)
{
    if(is_x_chr) // no conversion needed
        return NumericMatrix(0,0);

    NumericMatrix result(4,2);
    result(0,0) = 1.0;
    result(1,0) = result(1,1) = 0.5;
    result(2,0) = result(2,1) = 0.5;
    result(3,1) = 1.0;

    return result;
}

// check that sex conforms to expectation
const bool F2PK::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n = is_female.size();
    if(!any_x_chr) { // all autosomes
        if(n > 0) {
            result = true; // don't call this an error
            //REprintf("is_female included but not needed without X chromosome\n");
        }
    }
    else { // X chr included
        if(n == 0) {
            result = false;
            //REprintf("is_female not provided, but needed to handle X chromosome\n");
        }
        else {
            unsigned int n_missing = 0;
            for(unsigned int i=0; i<n; i++)
                if(is_female[i] == NA_LOGICAL) ++n_missing;
            if(n_missing > 0) {
                result = false;
                //REprintf("%d missing is_female values; is_female should not be missing.\n", n_missing);
            }
        }
    }
    return result;
}

// check that cross_info conforms to expectation
const bool F2PK::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n_row = cross_info.rows();
    const unsigned int n_col = cross_info.cols();

    if(!any_x_chr) { // all autosomes
        if(n_col > 0) {
            result = true; // don't call this an error
            //REprintf("cross_info included but not needed without X chromosome\n");
        }
    }
    else { // X chr included
        if(n_col == 0) {
            result = false;
            //REprintf("cross_info not provided, but needed to handle X chromosome\n");
        }
        else if(n_col > 1) {
            result = false;
            //REprintf("cross_info has %d columns, but should have just 1\n", n_col);
        }
        else {
            unsigned int n_missing = 0;
            for(unsigned int i=0; i<n_row; i++)
                if(cross_info[i] == NA_INTEGER) ++n_missing;
            if(n_missing > 0) {
                result = false;
                //REprintf("%d missing cross_info values; cross_info should not be missing.\n", n_missing);
            }

            unsigned int n_invalid = 0;
            for(unsigned int i=0; i<n_row; i++)
                if(cross_info[i] != NA_INTEGER && cross_info[i] != 0 && cross_info[i] != 1) ++n_invalid;
            if(n_invalid > 0) {
                result = false;
                //REprintf("%d invalid cross_info values; cross_info should be 0 or 1.\n", n_invalid);
            }
        }
    }
    return result;
}

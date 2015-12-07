// phase-known AIL QTLCross class (for HMM, in particular est.map)

#include "cross_ailpk.h"
#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "r_message.h"

enum gen {AA=1, AB=2, BA=3, BB=4,
          A=1, H=2, B=3, notB=4, notA=5,
          AY=5, BY=6};

const bool AILPK::check_geno(const int gen, const bool is_observed_value,
                            const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    // allow any value 0-5 or observed
    if(is_observed_value) {
        if(gen==0 || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        else return false;
    }

    if(is_x_chr && !is_female) {
        if(gen==AY || gen==BY) return true;
    }
    else { // autosome or female X
        if(gen==AA || gen==AB || gen==BA || gen==BB) return true;
    }

    return false; // invalid
}


const double AILPK::init(const int true_gen,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const int n_gen = cross_info[0];
    const int dir   = cross_info[1];

    if(is_x_chr) {
        if(dir==2) { // balanced case
            if(is_female) {
                return log(0.25);
            }
            else { // male
                return log(0.5);
            }
        }
        else { // AxB or BxA
            // frequency of A in AxB in males is (2/3) + (1/3)*(-1/2)^(s-1)
            double logm, logf, log1mm, log1mf;
            if(n_gen % 2 == 1) { // s is odd
                logf = log(2.0/3.0) + log1p( -exp( -((double)(n_gen+1) * log(2.0)) ) );
                logm = log(2.0/3.0) + log1pexp( -((double)(n_gen) * log(2.0)) );
            }
            else {
                logf = log(2.0/3.0) + log1pexp( -((double)(n_gen+1) * log(2.0)) );
                logm = log(2.0/3.0) + log1p( -exp( -((double)(n_gen) * log(2.0)) ) );
            }

            if(dir == 0) { // AxB
                log1mf = log1p(-exp(logf));
                log1mm = log1p(-exp(logm));
            }
            else { // BxA, take 1-p
                log1mf = logf;
                log1mm = logm;
                logf = log1p(-exp(logf));
                logm = log1p(-exp(logm));
            }

            if(is_female) {
                if(true_gen == AA) return 2.0*logf;
                else if(true_gen == AB || true_gen==BA) return logf + log1mf;
                else return 2.0*log1mf;
            }
            else { // male
                if(true_gen==AY) return logm;
                else return log1mm;
            }
        }
    }
    else {
        return log(0.25);
    }
}

const double AILPK::emit(const int obs_gen, const int true_gen, const double error_prob,
                        const IntegerVector& founder_geno, const bool is_x_chr,
                        const bool is_female, const IntegerVector& cross_info)
{

    #ifndef NDEBUG
    if(!check_geno(true_gen, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(obs_gen==0 || !check_geno(obs_gen, true, is_x_chr, is_female, cross_info))
        return 0.0; // missing or invalid

    if(is_female || !is_x_chr) { // female X chromosome just like autosome
        if(true_gen==AA) {
            if(obs_gen==A) return log(1.0-error_prob);
            if(obs_gen==H || obs_gen==B) return log(error_prob/2.0);
            if(obs_gen==notB) return log(1.0-error_prob/2.0);
            if(obs_gen==notA) return log(error_prob);
        }
        if(true_gen==AB || true_gen==BA) {
            if(obs_gen==H) return log(1.0-error_prob);
            if(obs_gen==A || obs_gen==B) return log(error_prob/2.0);
            if(obs_gen==notB || obs_gen==notA) return log(1.0-error_prob/2.0);
        }
        if(true_gen==BB) {
            if(obs_gen==B) return log(1.0-error_prob);
            if(obs_gen==H || obs_gen==A) return log(error_prob/2.0);
            if(obs_gen==notA) return log(1.0-error_prob/2.0);
            if(obs_gen==notB) return log(error_prob);
        }
    }
    else { // males
        if(true_gen==AY) {
            if(obs_gen==A || obs_gen==notB) return log(1.0-error_prob);
            if(obs_gen==B || obs_gen==notA) return log(error_prob);
            return 0.0; // treat everything else as missing
        }
        if(true_gen==BY) {
            if(obs_gen==B || obs_gen==notA) return log(1.0-error_prob);
            if(obs_gen==A || obs_gen==notB) return log(error_prob);
            return 0.0; // treat everything else as missing
        }
    }

    return NA_REAL; // shouldn't get here
}


const double AILPK::step(const int gen_left, const int gen_right, const double rec_frac,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const int n_gen = cross_info[0];
    const int dir = cross_info[1];

    if(is_x_chr) {
        if(dir == 2) { // balanced case
            double z = sqrt((1.0-rec_frac)*(9.0-rec_frac));
            double w = (1.0 - rec_frac + z)/4.0;
            double y = (1.0 - rec_frac - z)/4.0; // make positive
            double Rm, Rf;
            double k = (double)(n_gen - 2);
            double wk = pow(w,k);
            double yk = pow(y,k);
            Rm = 1 - 0.25*(2.0 + (1.0-2.0*rec_frac) * (wk + yk) +
                           (3.0 - 5.0*rec_frac + 2.0*rec_frac*rec_frac)/z * (wk - yk));
            Rf = 1 - 0.25*(2.0 + (1.0-2.0*rec_frac) * (wk + yk) +
                           (3.0 - 6.0*rec_frac + rec_frac*rec_frac)/z * (wk - yk));

            if(is_female) {
                switch(gen_left) {
                case AA:
                    switch(gen_right) {
                    case AA: return 2.0*log1p(-Rf);
                    case AB: case BA: return log1p(-Rf)+log(Rf);
                    case BB: return 2.0*log(Rf);
                    }
                case AB:
                    switch(gen_right) {
                    case AA: case BB: return log(Rf)+log1p(-Rf);
                    case AB: return 2.0*log(1.0-Rf);
                    case BA: return 2.0*log(Rf);
                    }
                case BA:
                    switch(gen_right) {
                    case AA: case BB: return log(Rf)+log1p(-Rf);
                    case BA: return 2.0*log(1.0-Rf);
                    case AB: return 2.0*log(Rf);
                    }
                case BB:
                    switch(gen_right) {
                    case AA: return 2.0*log(Rf);
                    case AB: case BA: return log1p(-Rf)+log(Rf);
                    case BB: return 2.0*log1p(-Rf);
                    }
                }
            }
            else { // male
                if(gen_left == gen_right) return log1p(-Rm);
                return log(Rm);
            }
        }
        else { // 0 = AxB; 1 = BxA
            // calculate frequency of AA haplotype in males and females
            double m11prev = 1.0;
            double f11prev = 0.5;
            double m11=1.0, f11=0.5, qp, qpp;
            for(int i=2; i<=n_gen; i++) {
                qpp = (2.0/3.0) + (1.0/3.0)*pow(-0.5, (double)(i-3));
                qp  = (2.0/3.0) + (1.0/3.0)*pow(-0.5, (double)(i-2));
                m11 = (1.0-rec_frac)*f11prev + rec_frac*qp*qpp;
                f11 = m11prev/2.0 + (1.0-rec_frac)/2.0*f11prev + (rec_frac/2.0)*qp*qpp;
                m11prev = m11;
                f11prev = f11;
            }

            if(is_female) {
                // allele frequencies at this generation
                double qf = (2.0/3.0) + (1.0/3.0)*pow(-0.5, (double)n_gen);
                // conditional probabilities along random haplotypes
                double f1to1 = f11/qf;
                double f1to2 = 1.0 - f1to1;
                double f2to1 = (qf - f11)/(1.0-qf);
                double f2to2 = 1.0 - f2to1;

                if(dir == 0) { // AxB
                    switch(gen_left) {
                    case AA:
                        switch(gen_right) {
                        case AA: return 2.0*log(f1to1);
                        case AB: case BA: return log(f1to1)+log(f1to2);
                        case BB: return 2.0*log(f1to2);
                        }
                    case AB:
                        switch(gen_right) {
                        case AA: return log(f1to1)+log(f2to1);
                        case AB: return log(f1to1)+log(f2to2);
                        case BA: return log(f1to2)+log(f2to1);
                        case BB: return log(f1to2)+log(f2to2);
                        }
                    case BA:
                        switch(gen_right) {
                        case AA: return log(f1to1)+log(f2to1);
                        case BA: return log(f1to1)+log(f2to2);
                        case AB: return log(f1to2)+log(f2to1);
                        case BB: return log(f1to2)+log(f2to2);
                        }
                    case BB:
                        switch(gen_right) {
                        case AA: return 2.0*log(f2to1);
                        case AB: case BA: return log(f2to2)+log(f2to1);
                        case BB: return 2.0*log(f2to2);
                        }
                    }
                }
                else { // BxA
                    switch(gen_left) {
                    case BB:
                        switch(gen_right) {
                        case BB: return 2.0*log(f1to1);
                        case AB: case BA: return log(f1to1)+log(f1to2);
                        case AA: return 2.0*log(f1to2);
                        }
                    case AB:
                        switch(gen_right) {
                        case BB: return log(f1to1)+log(f2to1);
                        case AB: return log(f1to1)+log(f2to2);
                        case BA: return log(f1to2)+log(f2to1);
                        case AA: return log(f1to2)+log(f2to2);
                        }
                    case BA:
                        switch(gen_right) {
                        case BB: return log(f1to1)+log(f2to1);
                        case BA: return log(f1to1)+log(f2to2);
                        case AB: return log(f1to2)+log(f2to1);
                        case AA: return log(f1to2)+log(f2to2);
                        }
                    case AA:
                        switch(gen_right) {
                        case BB: return 2.0*log(f2to1);
                        case AB: case BA: return log(f2to2)+log(f2to1);
                        case AA: return 2.0*log(f2to2);
                        }
                    }
                }
            }
            else {
                // allele frequencies at this generation
                double qm = (2.0/3.0) + (1.0/3.0)*pow(-0.5, (double)(n_gen-1));
                // conditional probabilities along random haplotypes
                double m1to1 = m11/qm;
                double m1to2 = 1.0 - m1to1;
                double m2to1 = (qm - m11)/(1.0-qm);
                double m2to2 = 1.0 - m2to1;

                if(dir == 0) { // AxB
                    switch(gen_left) {
                    case AY:
                        switch(gen_right) {
                        case AY: return log(m1to1);
                        case BY: return log(m1to2);
                        }
                    case BY:
                        switch(gen_right) {
                        case AY: return log(m2to1);
                        case BY: return log(m2to2);
                        }
                    }
                }
                else { // BxA; change A->B and B->A
                    switch(gen_left) {
                    case AY:
                        switch(gen_right) {
                        case AY: return log(m2to2);
                        case BY: return log(m2to1);
                        }
                    case BY:
                        switch(gen_right) {
                        case AY: return log(m1to2);
                        case BY: return log(m1to1);
                        }
                    }
                }


            }
        }
    }
    else { // autosome
        // cross_info[0] is the number of generations
        // R = 0.5*[ 1 - (1-2r)*(1-r)^(s-2) ]
        double tmp = (1-2.0*rec_frac)*pow(1.0-rec_frac, (double)(n_gen-2));
        double logR = -log(2.0) + log1p(-tmp);
        double log1mR = -log(2.0) + log1p(tmp);

        switch(gen_left) {
        case AA:
            switch(gen_right) {
            case AA: return 2.0*log1mR;
            case AB: case BA: return log1mR+logR;
            case BB: return 2.0*logR;
            }
        case AB:
            switch(gen_right) {
            case AA: case BB: return logR+log1mR;
            case AB: return 2.0*log1mR;
            case BA: return 2.0*logR;
            }
        case BA:
            switch(gen_right) {
            case AA: case BB: return logR+log1mR;
            case BA: return 2.0*log1mR;
            case AB: return 2.0*logR;
            }
        case BB:
            switch(gen_right) {
            case AA: return 2.0*logR;
            case AB: case BA: return log1mR+logR;
            case BB: return 2.0*log1mR;
            }
        }
    }

    return NA_REAL; // shouldn't get here
}

const IntegerVector AILPK::possible_gen(const bool is_x_chr, const bool is_female,
                                       const IntegerVector& cross_info)
{
    if(is_x_chr && !is_female) { // male X chromosome
        IntegerVector result = IntegerVector::create(AY,BY);
        return result;
    }
    else { // autosome
        IntegerVector result = IntegerVector::create(AA,AB,BA,BB);
        return result;
    }
}

const int AILPK::ngen(const bool is_x_chr)
{
    if(is_x_chr) return 6;
    return 4;
}

const double AILPK::nrec(const int gen_left, const int gen_right,
                        const bool is_x_chr, const bool is_female,
                        const IntegerVector& cross_info)
{
    #ifndef NDEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_x_chr && !is_female) { // male X chr
        if(gen_left == gen_right) return 0.0;
        else return 1.0;
    }
    else { // autosome or female X chr
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

const double AILPK::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                const IntegerMatrix& cross_info, const int n_gen)
{
    Rcpp::stop("est_map not yet available for AIL");

    return NA_REAL;
}


const NumericMatrix AILPK::geno2allele_matrix(const bool is_x_chr)
{
    if(is_x_chr) {
        NumericMatrix result(6,2);

        result(0,0) = 1.0;               // AA female
        result(1,0) = result(1,1) = 0.5; // AB female
        result(2,0) = result(2,1) = 0.5; // BA female
        result(3,1) = 1.0;               // BB female

        result(4,0) = 1.0; // AY male
        result(5,1) = 1.0; // BY male

        return result;
    }
    else {
        NumericMatrix result(4,2);
        result(0,0) = 1.0;
        result(1,0) = result(1,1) = 0.5;
        result(2,0) = result(2,1) = 0.5;
        result(3,1) = 1.0;

        return result;
    }
}

// check that sex conforms to expectation
const bool AILPK::check_is_female_vector(const LogicalVector& is_female, const bool any_x_chr)
{
    bool result = true;
    const unsigned int n = is_female.size();
    if(!any_x_chr) { // all autosomes
        if(n > 0) {
            // not needed here, but don't call this an error
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
const bool AILPK::check_crossinfo(const IntegerMatrix& cross_info, const bool any_x_chr)
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

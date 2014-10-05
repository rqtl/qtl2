// phase-known intercross HMM functions (for est.map)

#include <math.h>
#include <Rcpp.h>
#include "cross.h"

enum gen {GENO_NA=0, AA=1, AB=2, BA=3, BB=4,
          A=1, H=2, B=3, notB=4, notA=5,
          AY=1, BY=4};

bool F2PK::check_geno(int gen, bool is_observed_value,
                      bool is_X_chr, bool is_female, IntegerVector cross_info)
{
    // allow any value 0-5 or observed
    if(is_observed_value) {
        if(gen==GENO_NA || gen==A || gen==H || gen==B ||
           gen==notA || gen==notB) return true;
        throw std::range_error("Invalid observed genotype");
    }

    if(is_X_chr) {
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

    throw std::range_error("Invalid genotype");
    return false; // can't get here
}


double F2PK::init(int true_gen,
                  bool is_X_chr, bool is_female,
                  IntegerVector cross_info)
{
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    if(is_X_chr) return log(0.5);
    else return log(0.25);
}

double F2PK::emit(int obs_gen, int true_gen, double error_prob,
                  bool is_X_chr, bool is_female,
                  IntegerVector cross_info)
{
    check_geno(obs_gen, true, is_X_chr, is_female, cross_info);
    check_geno(true_gen, false, is_X_chr, is_female, cross_info);

    if(obs_gen==GENO_NA) return 0.0; // log(1.0)

    if(is_X_chr) {
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

    throw std::range_error("invalid obs_gen or true_gen");
    return NA_REAL; // can't get here
}


double F2PK::step(int gen_left, int gen_right, double rec_frac,
                  bool is_X_chr, bool is_female,
                  IntegerVector cross_info)
{
    check_geno(gen_left, false, is_X_chr, is_female, cross_info);
    check_geno(gen_right, false, is_X_chr, is_female, cross_info);

    if(is_X_chr) {
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

    throw std::range_error("invalid genotypes.");
    return NA_REAL; // can't get here
}

IntegerVector F2PK::possible_gen(bool is_X_chr, bool is_female,
                                 IntegerVector cross_info)
{
    if(is_X_chr) {
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
            int vals[] = {AA,BB};
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

int F2PK::ngen(bool is_X_chr)
{
    return 4;
}

double F2PK::nrec(int gen_left, int gen_right,
                bool is_X_chr, bool is_female,
                IntegerVector cross_info)
{
    check_geno(gen_left, false, is_X_chr, is_female, cross_info);
    check_geno(gen_right, false, is_X_chr, is_female, cross_info);

    if(is_X_chr) {
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

    throw std::range_error("invalid genotypes.");
    return NA_REAL; // can't get here
}

double F2PK::est_rec_frac(NumericMatrix gamma, bool is_X_chr)
{
    int n_gen = gamma.rows();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0;
    for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

    if(is_X_chr) {
        double diagsum = 0.0;
        for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

        return(1.0 - diagsum/denom);
    }

    IntegerVector empty(0);

    double numerator = 0.0;
    for(int il=0; il<n_gen; il++)
        for(int ir=0; ir<n_gen; ir++)
            numerator += gamma(il,ir)*nrec(il+1, ir+1, false, false, empty);

    return numerator/denom;
}

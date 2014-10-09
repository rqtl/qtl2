// RI sib-mating HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_risib.h"

enum gen {AA=1, BB=2};

double RISIB::init(int true_gen,
                   bool is_X_chr, bool ignored, IntegerVector cross_info)
{
    #ifdef DEBUG
    if(!check_geno(true_gen, false, is_X_chr, ignored, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const bool forward_direction = (cross_info[0] == 0); // AA female x BB male

    if(is_X_chr) {
        if(forward_direction) {
            if(true_gen == AA) return log(2.0)-log(3.0);
            if(true_gen == BB) return -log(3.0);
        }
        else {
            if(true_gen == BB) return log(2.0)-log(3.0);
            if(true_gen == AA) return -log(3.0);
        }
    }
    else { // autosome
        return -log(2.0);
    }

    return NA_REAL; // can't get here
}

double RISIB::step(int gen_left, int gen_right, double rec_frac,
                   bool is_X_chr, bool ignored, IntegerVector cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_X_chr, ignored, cross_info) ||
       !check_geno(gen_right, false, is_X_chr, ignored, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    if(is_X_chr)  {
        const bool forward_direction = (cross_info[0] == 0); // AA female x BB male

        if(forward_direction) {
            if(gen_left == AA) {
                if(gen_right == AA)
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == BB)
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }

            if(gen_left == BB) {
                if(gen_right == BB)
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == AA)
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
        else {
            if(gen_left == AA) {
                if(gen_right == AA)
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == BB)
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }

            if(gen_left == BB) {
                if(gen_right == BB)
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == AA)
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
    }
    else { // autosome
        const double R = 4.0*rec_frac/(1+6.0*rec_frac);

        if(gen_left == gen_right) return log(1.0-R);
        else return log(R);
    }

    return NA_REAL; // can't get here
}

double RISIB::est_rec_frac(NumericMatrix gamma, bool is_X_chr)
{
    int n_gen = gamma.rows();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0;
    for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

    double diagsum = 0.0;
    for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

    double R = 1.0 - diagsum/denom;

    return R/(4.0-6.0*R);
}

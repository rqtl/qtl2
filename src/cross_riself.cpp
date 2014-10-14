// RI selfing HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_riself.h"

enum gen {AA=1, BB=2};

const double RISELF::step(const int gen_left, const int gen_right, const double rec_frac,
                          const bool is_x_chr, const bool is_female, const IntegerVector& cross_info)
{
    #ifdef DEBUG
    if(!check_geno(gen_left, false, is_x_chr, is_female, cross_info) ||
       !check_geno(gen_right, false, is_x_chr, is_female, cross_info))
        throw std::range_error("genotype value not allowed");
    #endif

    const double R = 2.0*rec_frac/(1+2.0*rec_frac);

    if(gen_left == gen_right) return log(1.0-R);
    else return log(R);
}

const double RISELF::est_rec_frac(const NumericVector& gamma, const bool is_x_chr,
                                  const IntegerMatrix& cross_info, const int n_gen)
{
    int n_ind = cross_info.cols();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0, diagsum = 0.0;
    for(int ind=0, offset=0; ind<n_ind; ind++, offset += n_gen_sq) {
        for(int i=0; i<n_gen_sq; i++) denom += gamma[offset+i];
        for(int i=0; i<n_gen; i++)    diagsum += gamma[offset+i*n_gen+i];
    }

    double R = 1.0 - diagsum/denom;

    return 0.5*R/(1-R);
}

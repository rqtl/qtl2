// RI selfing HMM functions

#include <math.h>
#include <Rcpp.h>
#include "cross.h"
#include "cross_riself.h"

enum gen {AA=1, BB=2};

double RISELF::step(int gen_left, int gen_right, double rec_frac,
                    bool ignored1, bool ignored2, IntegerVector ignored3)
{
    check_geno(gen_left, false, ignored1, ignored2, ignored3);
    check_geno(gen_right, false, ignored1, ignored2, ignored3);

    const double R = 2.0*rec_frac/(1+2.0*rec_frac);

    if(gen_left == gen_right) return log(1.0-R);
    else return log(R);
}

double RISELF::est_rec_frac(NumericMatrix gamma, bool is_X_chr)
{
    int n_gen = gamma.rows();
    int n_gen_sq = n_gen*n_gen;

    double denom = 0.0;
    for(int i=0; i<n_gen_sq; i++) denom += gamma[i];

    double diagsum = 0.0;
    for(int i=0; i<n_gen; i++) diagsum += gamma(i,i);

    double R = 1.0 - diagsum/denom;

    return 0.5*R/(1-R);
}

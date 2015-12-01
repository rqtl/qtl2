// genome scan by Haley-Knott regression

#include "scan_hk.h"
#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"

// Scan a single chromosome with no additive covariates (not even intercept)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_nocovar(NumericVector genoprobs, NumericMatrix pheno,
                                     const double tol=1e-12)
{
    const unsigned int n_ind = pheno.rows();
    const unsigned int n_phe = pheno.cols();
    Dimension d = genoprobs.attr("dim");
    const unsigned int n_pos = d[2];
    const unsigned int n_gen = d[1];
    const unsigned int x_size = n_ind * n_gen;
    // check that d[0] == n_ind;

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen);

    for(unsigned int i=0, offset=0; i<n_pos; i++, offset += x_size) {
        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // calc rss and paste into ith column of result
        result(_,i) = calc_rss_linreg(X, pheno, tol);
    }

    return result;
}

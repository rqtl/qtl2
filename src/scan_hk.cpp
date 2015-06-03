// genome scan by Haley-Knott regression

#include <Rcpp.h>

using namespace Rcpp;

#include "linreg.h"
#include "scan_hk.h"

// Scan a single chromosome with no additive covariates (not even intercept)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
//
// [[Rcpp::export]]
NumericMatrix scan_hk_onechr_nocovar(NumericVector genoprobs, NumericMatrix pheno)
{
    int n_ind = pheno.rows();
    int n_phe = pheno.cols();
    Dimension d = genoprobs.attr("dim");
    int n_pos = d[2];
    int n_gen = d[1];
    int x_size = n_ind * n_gen;
    // check that d[0] == n_ind;

    NumericMatrix result(n_phe, n_pos);
    NumericMatrix X(n_ind, n_gen);

    for(int i=0, offset=0; i<n_pos; i++, offset += x_size) {
        // copy genoprobs for pos i into a matrix
        std::copy(genoprobs.begin() + offset, genoprobs.begin() + offset + x_size, X.begin());

        // calc rss and paste into ith column of result
        result(_,i) = calc_rss_linreg(X, pheno);
    }

    return result;
}

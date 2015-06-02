// genome scan by Haley-Knott regression

// Scan a single chromosome with no additive covariates (not even intercept)
//
// genoprobs = 3d array of genotype probabilities (individuals x genotypes x positions)
// pheno     = matrix of numeric phenotypes (individuals x phenotypes)
//             (no missing data allowed)
//
// output    = matrix of residual sums of squares (RSS) (phenotypes x positions)
NumericMatrix scan_hk_onechr_nocovar(const NumericVector& genoprobs, const NumericMatrix& pheno);

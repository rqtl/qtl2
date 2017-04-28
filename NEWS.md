## qtl2scan 0.5-12 (2017-04-28)

### New features

- Implemented model="binary" (for phenotypes with values 0/1) in scan1,
  scan1coef, and fit1. (Not available with kinship correction.)


## qtl2scan 0.5-11 (2017-04-27)

### Minor fixes

- Reduce repeated code dealing with kinship matrices.

### Bug fixes

- Fix bug regarding treatment of pre-decomposed kinship matrix in
  `scan1`.


## qtl2scan 0.5-9 (2017-04-19)

### Minor changes

- Revised `subset_scan1` (and the internal function `subset_kinship`)
  to use the same options for subsetting by chromosome as the
  functions in [R/qtl2geno](https://github.com/rqtl/qtl2geno), most
  importantly use of "negative" chromosome indexes, like `"-X"`.


## qtl2scan 0.5-8 (2017-04-03)

### New features

- Added `fit1()` to fit a single-QTL model at a single fixed position
  and return the LOD score, estimated coefficients, individual
  contributions to the LOD score, and (if `se=TRUE`) standard errors.

- In `scan1coef()` and `scan1blup()`, added an argument `nullcovar`
  for covariates to include only under the null hypothesis (of no
  QTL). This is only used in the case that `kinship` is provided but
  `hsq` is not, as these may be needed for the X chromosome to get the
  estimated residual heritability.


## qtl2scan 0.5-7 (2017-03-21)

### New features

- Added `dim.calc_genoprob` and `dimnames.calcgenoprob`, from
  [Brian Yandell](https://github.com/byandell), for use with
  [qtl2feather](https://github.com/byandell/qtl2feather), which uses
  [feather](https://github.com/wesm/feather) to store
  genotype probabilities in a file (to save memory).

- In precess of revising various functions to use qtl2feather,
  particularly in grabbing dimnames (with the above functions), but
  also to avoid `seq(along=genoprobs)` and instead use
  `seq_len(length(genoprobs))`.


## qtl2scan 0.5-6 (2017-03-17)

### New features

- Added `scan1perm` to perform a permutation test to establish
  genome-wide significance in a single-QTL genome scan by `scan1`.

- Also added functions `rbind.scan1perm` and `cbind.scan1perm` for
  combining `scan1perm` results, and `summary.scan1perm` to obtain
  significance thresholds.

### Bug fixes

- Fixed a bug in `scan1`. This would only show up if you were using a
  kinship matrix and scanning the X chromosome on its own with special
  X chr covariates (`Xcovar`). (Accidentally was acting as if it were
  an autosome and so ignoring `Xcovar`.)


## qtl2scan 0.5-5 (2017-03-13)

### New features

- Refactored to simplify the main data structures for `scan1`,
  `scan1coef`, and `scan1blup` output, and to deal with the
  refactoring of data structures in [qtl2geno](https://github.com/rqtl/qtl2geno).
  Functions like `max_scan1`, `find_peaks`, `lod_int`, and `bayes_int`
  now need you to provide a map. Similarly, to subset `scan1` results
  by chromosome, you need to provide a map to the subsetting function.

- Pulled the `"snpinfo"` attribute out of the `scan1` object. Now you
  need to use `index_snps` to identify groups of equivalent SNPs prior
  to running `genoprob_to_snpprob`. `index_snps` adds some new columns
  to the `snpinfo` data frame, which are then needed by `top_snps`.

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

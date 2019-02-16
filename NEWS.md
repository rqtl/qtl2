## qtl2 0.19-2 (2019-02-16)

### Bug fixes

- Fix bug in `zip_datafiles()` when the files are in a subdirectory.
  (See [Issue #102](https://github.com/rqtl/qtl2/issues/102).)


## qtl2 0.18 (2019-02-08)

### New features

- Added `plot_lodpeaks()` for scatterplot of LOD score vs position for
  inferred QTL from `find_peaks()` output.

- Added new cross types `"genril"` and `"genail"`, implemented to
  handle any number of founders; include the number of founders in the
  cross type, for example `"genril38"` or `"genail38"`. The cross
  information has length 1 + number of founders, with first column
  being the number of generations and the remaining columns being
  non-negative integers that indicate the relative frequencies of the
  founders in the initial population (these will be scaled to sum to
  1). `"genril"` assumes the progeny are inbred lines (recombinant
  inbred lines, RIL), while `"genail"` assumes the progeny have two
  random chromosomes (advanced intercross lines, AIL).

- The internal function `batch_vec()` now made user-accessible, and
  takes an additional argument `n_cores`. This splits a vector into
  batches for use in parallel calculations.

- The internal function `cbind_expand()` now made user-accessible.
  It's for combining matrices using row names to align the rows and
  expanding with missing values if there are rows in some matrices but
  not others.

- In `plot_peaks()`, added `lod_labels` argument. If TRUE, include LOD
  scores as text labels in the figure.

- Added function `calc_het()` for calculating estimated
  heterozygosities, by individual or by marker, from genotype
  probabilities derived by `calc_genoprob()`.

### Minor changes

- Small corrections to documentation.

- Revise some tests due to change in Recla and DOex datasets at
  <https://github.com/rqtl/qtl2data>

- Add tests of decomposed kinship matrix (from `decomp_kinship()`)
  with `scan1()`.

- `rbind_scan1()` and `cbind_scan1()` no longer give error if inputs
  don't all have matching attributes.

- Change default gap between chromosomes in `plot_scan1()` (and
  related) to be 1% of the total genome length.

### Bug fixes

- Fixed bug in `subset_kinship()` that prevented `scan1()` from
  working with decomposed "loco" kinship matrices.

- Fixed descriptions in help files for `cbind.calc_genoprob()` and
  `rbind.calc_genoprob()`, for column- and row-binding genotype
  probabilities objects (as output by `calc_genoprob()`. `cbind()` is
  for the same set of individuals but different chromosomes. `rbind()`
  is for the same set of markers and genotypes but different
  individuals. Made similar corrections for the related functions for
  `sim_geno()` and `viterbi()` output.


## qtl2 0.16 (2018-07-23)

### New features

- Added `pull_genoprobint()` for pulling out the genotype
  probabilities for a given genomic interval. Useful, for example, to
  apply `scan1blup()` over a defined interval rather than an entire
  chromosome.

- `scan1()`, `scan1perm()`, `scan1coef()`, `fit1()`, and `scan1snps()`
  can now use weights when `kinship` is provided, for example for the
  case of the analysis of recombinant inbred line (RIL) phenotype
  means with differing numbers of individuals per line. The residual
  variance matrix is like $v[h^2 K + (1-h^2)D]$ where D is diagonal
  {1/w} for weights w.

- Add `weights` argument to `est_herit()`.

- Added `add_threshold()` for adding significance thresholds to a
  genome scan plot.

- Added `predict_snpgeno()` for predicting SNP genotypes in a
  multiparent populations, from inferred genotypes plus the founder
  strains' SNP alleles.

- In `genoprob_to_snpprob()`, the `snpinfo` argument can now be a
  cross object (for a multiparent population with founder genotypes),
  in which case the SNP information for all SNPs with complete founder
  genotype data is calculated and used.

- `max_scan1()` with `lodcolumn=NULL` returns the maximum for all
  lod score columns. If `map` is included, the return value is in the
  form returned by `find_peaks()`, namely with `lodindex` and
  `lodcolumn` arguments added at the beginning.

- Added `replace_ids()` for replacing individual IDs in an object.
  S3 method for `"cross2"` objects and output of `calc_genoprob()`,
  `viterbi()`, `maxmarg()`, and `sim_geno()`.

- Added `clean_scan1()` plus generic function `clean()` that works
  with both this and with `clean_genoprob()`. `clean_scan1()` replaces
  negative values with `NA` and removes rows that have all `NAs`.

### Minor changes

- More informative error message in `est_herit()`, `scan1()`, etc.,
  when covariates and other data are not numeric.

- Fixed `pull_genoprobpos()` so it will work with
  [qtl2feather](https://github.com/byandell/qtl2feather)
  (and [qtl2fst](https://github.com/rqtl/qtl2fst)).

- In `plot_genes()`, if `xlim` is provided as an argument, subset the
  genes to those that will actually appear in the plotting region.

- Revise `find_marker()` so that the input `map` can also be a "snp
  info" table (with columns `"snp_id"`, `"chr"` and `"pos"`).

- Added `find_index_snp()` for identifying the index SNP that
  corresponds to a particular SNP in a snp info table that's been
  indexed with `index_snps()`.

- Add `overwrite` argument (default `FALSE`) to `zip_datafiles()`,
  similar to that for `write_control_file()`.

- `plot_snpasso()` now takes an argument `chr`.

- `max_scan1()` no longer gives a warning if `map` is not provided.

- `insert_pseudomarkers()` will now accept `pseudomarker_map` that
   includes only a portion of the chromosomes.

- In `fit1()`, replaced `tol` and `maxit` and added `...` which takes
  these plus a few additional hidden control parameters.

### Bug fixes

- Fix a bug in `index_snps()`; messed up results when `start` and
  `end` outside the range of the map.

- Fix a bug in `scan1snps()` regarding use of `chr` argument: need to
  force to be unique character strings, and avoid unnecessary warning
  about `start` and `end`.

- Fix a bug in `scan1snps()` where it didn't check that the `genoprobs`
  and `map` conform.

- Revised underlying binary trait regression function to avoid some of
  the tendency towards NAs.


## qtl2 0.14 (2018-03-09)

### New features

- Added function `clean_genoprob()` which cleans genotype
  probabilities by setting small values to 0 and, for genotype columns
  where the maximum value is not large, setting all values to 0. This
  is intended to help with the problem of unstable estimates of
  genotype effects in `scan1coef()` and `fit1()` when there's a
  genotype that is largely absent.

- Added function `compare_maps()` for comparing marker order between
  two marker maps.

- Revised the order of arguments in `reduce_markers()` to match
  `pick_marker_subset()`, because I like the latter better. Removed
  the function `pick_marker_subset()` because it's identical to
  `reduce_markers()`. (Seriously, I implemented the same thing twice.)

### Minor changes

- `plot_coef()` now uses a named `lodcolumn` argument, if provided, to
  subset `scan1_output`, if that's provided.

- In the documentation for `scan1coef()`, `scan1blup()`, and `fit1()`,
  revised the suggested contrasts for getting additive and dominance
  effects in an intercross.


### Bug fixes

- In `plot_coef()` with `scan1_output` provided, `ylim_lod` was being
  ignored.


## qtl2 0.12 (2018-01-19)

### New features

- `find_peaks()` and `max_scan1()` can now take snpinfo tables (as
  produced by `index_snps()` and `scan1snps()`) in place of the map.

### Minor changes

- Sped up a bunch of the examples in the help files (mostly by
  subsetting the example datasets).

### Bug fixes

- Further embarassment: the bug fix in version 0.10 didn't fully fix
  the problem with `find_peaks()`.


## qtl2 0.10 (2018-01-09)

### Bug fixes

- Fixed embarassing bug in `find_peaks()`. (Stopped with error if no
  LOD scores were above the threshold.)


## qtl2 0.8 (2018-01-05)

- First formal release


## qtl2 0.7-6 (2017-12-12)

### New features

- The output of `fit1()` now includes fitted values.

- Added function `pull_genoprobpos()` for pulling out a specific
  position (by name or position) from a set of genotype probabilities.

- The `chr` column in the result of `find_peaks()` is now a factor.
  This makes it possible to sort by chromosome. Also added an argument
  `sort_by` for choosing how to sort the rows in the result (by
  column, genomic position, or LOD score).

- In `max_scan1()`, if `map` is *not* provided, rather than stopping
  with an error, we just issue a warning and return the genome-wide
  maximum LOD score.

- Revised `find_markerpos()` so it can take a map (as a list of
  vectors of marker positions) in place of a `"cross2"` object.

### Minor changes

- Added many more checks of the inputs to various functions.

### Bug fixes

- Fixed a bug in `plot.scan1()`, which failed to pass `lodcolumn` to
  `plot_snpasso()`.


## qtl2 0.7-1 (2017-11-27)

The previously separate packages qtl2geno, qtl2scan, qtl2plot, and
qtl2db have now been combined into one package.

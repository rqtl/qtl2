## qtl2 0.23-10 (2020-12-03)

### Major changes

- Revised treatment of X chromosome in "general AIL" cross type, to be
  the same as the autosomes but with 2/3 as many generations for
  recombination. This should provide a better approximation.

### Minor changes

- Revised `reduce_markers()` so that it can handle the case of _many_
  markers, by working with them in smaller batches.

- `fit1()` now returns both fitted values and residuals.

- Updated mouse gene database with 2020-09-07 data from
  [MGI](http://www.informatics.jax.org/downloads/mgigff3/archive/monthly/).

- Implemented Issue #184, to make `calc_het()` multi-core.

- Making the tests of the plot functions only run locally, as they are
  failing on R-devel and it seems like it could be a pain running
  these on CRAN.

### Bug fixes

- Fixed [Issue #181](https://github.com/rqtl/qtl2/issues/181), where
  `calc_het()` gave values > 1 when used with
  [R/qtl2fst](https://github.com/rqtl/qtl2fst)-based probabilities.
  Also fixed a similar bug in `calc_geno_freq()`.

- Fixed [Issue #172](https://github.com/rqtl/qtl2/issues/172), where
  `fit1()` gave incorrect fitted values when `kinship` is provided,
  because they weren't "rotated back".

- `fit1()` no longer provides individual LOD scores (`ind_lod`) when
  `kinship` is used, as with the linear mixed model, the LOD score is
  not simply the sum of individual contributions.

- Fixed [Issue #174](https://github.com/rqtl/qtl2/issues/174) re
  `genoprob_to_alleleprob()` for general AIL crosses. We had not
  implemented the `geno2allele_matrix()` function.

- Fix Issue #164, so `plot_pxg()` can handle a phenotype that is a
  single-column data frame.

- Fix Issue #135, so `plot_scan1()` can take vector input (which is
  then converted to a single-column matrix).


## qtl2 0.22-11 (2020-07-09)

### Bug fixes

- Small changes in tests to avoid failures on solaris


## qtl2 0.22-10 (2020-06-30)

### Bug fixes

- Fixed compilation error on Solaris on CRAN, due to a log(10) and
  sqrt(10). Solaris refuses to do log(int) or sqrt(int), I guess.

- Fixed some conflicting enum definitions in c++ code


## qtl2 0.22-8 (2020-06-22)

### Minor changes

- Added recognition of the R Core Team as a contributor, as the
  package includes a copy of code for Brent's method for univariate
  function optimization. Also added a Copyright field in the
  DESCRIPTION field, explaining the copyright of that code.

- Sped up some of the examples and tests. Tests no longer use more
  than 2 cores (even those that are only run locally).

- Added `\value{}` sections in the documentation for various
  functions. Added further explanation of `"viterbi"`, `"scan1"`,
  `"scan1perm"`, etc. S3 classes in the documentation.


## qtl2 0.22 (2020-05-21)

### Major changes

- Added some functions for diagnostics: `recode_snps()`,
  `calc_raw_het()`, `calc_raw_geno_freq()`, `calc_raw_maf()`, and
  `calc_raw_founder_maf()`.

- Added argument `blup` to `fit1()`, for getting BLUPs for a single
  fixed QTL position. At present, just gives estimates and
  coefficients by calling `scan1blup()` with a single position.

- `pull_genoprobpos()` can now take either a marker name (as before) or
  a set of map, chromosome, and position (from which it uses
  `find_marker()` to get the marker name).

- Added plot function for the results of `compare_geno()`. (Plots
  histogram of upper triangle.)

- Added functions `n_founders()` and `founders()` for getting the
  number of founders and the founder strain names for a cross2 object.

- `scan1()` now takes an optional `hsq` argument, so that the residual
  heritability may be specified rather than estimated.

### Minor changes

- `write_control_file()` now allows cross info codes with a cross info
  file (previously only allowed with a covariate). `read_cross2()`
  gives a warning if there are cross info conversion codes but more
  than one cross info column.

- Small fix in `read_cross2()` to allow multiple cross info covariates.

- Added a check that the founder genotypes have the same strain IDs on
  each chromosome.

- `convert2cross2()` now includes `alleles` component even if it
  wasn't present as an attribute.

- Added function `sdp2char()` for converting numeric SDP codes to
  character strings like `"ABC|DEFGH"`.

- Updated mouse gene database with 2019-08-12 data from
  [MGI](http://www.informatics.jax.org/downloads/mgigff3/archive/monthly/).

- `get_common_ids()` strips off names from output, just in case.

- Added internal functions `rqtl1_crosstype()` and `rqtl1_chrtype()`.

### Bug fixes

- Fixed typo in help for `scan1()` and related functions.

- `genoprob_to_snpprob()` was giving an error if you gave a cross2
  object in place of a snpinfo table and it had monomorphic markers.

- Fixed problem with weights in `scan1()` and related functions when
  their derived from `table()`. Make sure they're a plain numeric
  vector, not an array.

- Fixed `check_cross2()`: the check for invalid genotypes wasn't
  happening.

- Better error message for the case that there are no markers in
  common between map and genotypes.

- `extract_dim_from_header()`, used by `read_cross2()` and `read_csv()`,
  now just looks for the number part in the rest of the line.

- `maxlod()` now handles missing values (forcing `na.rm=TRUE`). If all
  values are missing it gives a warning and returns `-Inf`.
  [Fixes Issue #134.]

- In `max_scan1()`, treat the case that the input has no column names.
  [Fixes Issue #133.]

- `max_scan1()` was giving a messed up error message if `lodcolumn`
  was out of range. [Fixes Issue #132.]

- Revised the script `inst/scripts/create_ccvariants.R` to capture _all_
  of the consequences and genes for each SNP (rather than just the
  first), and fixing a bug that prevented capture of indels from
  chromosomes 6-X. Consequently, revised the example SQLite database
  `extdata/cc_variants_small.sqlite` and associated tests.


## qtl2 0.20 (2019-06-03)

### Major changes

- `scan1coef()` and `fit1()` now, by default, gives coefficient
  estimates for the QTL effects that sum to 0, with an additional
  coefficient being the intercept. This makes it more like DOQTL (and
  `scan1blup()`). The previous behavior can be obtained with the
  argument `zerosum=FALSE`.

- Added function `create_snpinfo()` for creating a SNP information table
  from a cross2 object, for use with `scan1snps()`.

### Minor changes

- Updated `extdata/mouse_genes_small.sqlite` using updated MGI
  annotations. Some of the field names have changed.

- In `check_cross2()`, added a test for alleles being a vector of
  character strings.

- Fixed some tests for R 3.6, due to change in random number generation.

- Use Markdown for function documentation, throughout

### Bug fixes

- In `genoprob_to_snpprob()` when a cross object is provided, make
  sure the genotype probabilities get subset to the cross markers.

- Fixed bug in `scan1snps()` re `keep_all_snps=FALSE`. It wasn't
  subsetting to the index SNPs properly. Added an internal function
  `reduce_to_index_snps()`.
  (See [Issue #89](https://github.com/rqtl/qtl2/issues/89).)

- Fixed bug in step probabilities for 4-, 8-, and 16-way RIL by selfing.

- Fixed bug in `zip_datafiles()` when the files are in a subdirectory.
  (See [Issue #102](https://github.com/rqtl/qtl2/issues/102).)

- Fixed bug in `plot_peaks()` for the case that the input `peaks`
  object does not contain QTL intervals.
  (See [Issue #107](https://github.com/rqtl/qtl2/issues/107).)

- Fixed inappropriate warning message for check of `cross_info` with
  cross type `risib8`.

- Fixed bugs in `guess_phase()` and `locate_xo()` where we needed an
  `any()` around a comparison of two vectors.


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

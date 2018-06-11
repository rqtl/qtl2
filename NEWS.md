## qtl2 0.15-15 (2018-06-11)

### New features

- Added `pull_genoprobint()` for pulling out the genotype
  probabilities for a given genomic interval. Useful, for example, to
  apply `scan1blup()` over a defined interval rather than an entire
  chromosome.

- `scan1()`, `scan1perm()`, and `scan1coef()` now can use weights when
  `kinship` is provided, for example for the case of the analysis of
  recombinant inbred line (RIL) phenotype means with differing numbers
  of individuals per line. The model is like $v[h^2 K + (1-h^2)D]$
  where D is diagonal {1/w} for weights w.

- Add `weights` argument to `est_herit()`.

- Added `add_threshold()` for adding significance thresholds to a
  genome scan plot.

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

### Bug fixes

- Fix a bug in `index_snps()`; messed up results when `start` and
  `end` outside the range of the map.


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

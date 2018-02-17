## qtl2 0.13-1 (2018-02-16)

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

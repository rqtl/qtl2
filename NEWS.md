## qtl2geno 0.5-24 (2017-06-12)

### Bug fixes

- Fix formatting problem in output of `summary.cross2()` within
  RStudio.


## qtl2geno 0.5-23 (2017-06-05)

### Minor changes

- Revised installation instructions.


## qtl2geno 0.5-22 (2017-05-24)

### New features

- Added function `scale_kinship()` which converts a kinship matrix (or
  a list of such, in the case of the "leave one chromosome out"
  method) to be like a correlation matrix.

- Removed the "normalize" argument from `calc_kinship()`, though left
  the internal function `normalize_kinship()` in place, for now.

### Minor changes

- `insert_pseudomarkers` now gives an error if the input `map` is
  `NULL`.


## qtl2geno 0.5-21 (2017-05-10)

### New features

- `count_xo` now works with the output of `sim_geno`. The result is a
  3d-array of counts of crossovers for each individual on each
  chromosome in each imputation.


## qtl2geno 0.5-20 (2017-05-04)

### New features

- Added `chisq_colpairs` which performs chi-square tests for
  independence on all pairs of columns of a matrix. It just calculates
  the statistics.


## qtl2geno 0.5-19 (2017-05-01)

### Minor changes

- Added an argument `save_rf` to `est_map()`; if `TRUE`, the estimated
  recombinations are saved as an attribute (`"rf"`) of the result.
  This can be useful for diagnostic purposes, for example when the
  estimated recombination fraction between markers is > 1/2.
  (After converting to genetic distance, rf>1/2 is indistinguishable
  from rf=1/2.)


## qtl2geno 0.5-18 (2017-04-19)

### New features

- New function `reduce_map_gaps` that reduces the length of any gaps
  in map. (Gaps greater than `min_gap` are reduced to `min_gap`.)

- `maxmarg` now picks at random among genotypes that jointly share the
  maximum probability. Previously, it picked the first among these.
  Added an argument `tol`; if two genotypes have probabilities that
  differ by no more than `tol`, they are treated as having the same
  probability.

- New function `calc_entropy` takes the results of `calc_genoprob` and
  calculates, for each individual at each genomic postion, the entropy
  of the genotype probability distribution, as a measure of missing
  information.

### Bug fixes

- Fix bug in `find_map_gaps` regarding the case that the output are
  empty.

- Fix bug in attempting to subsett `calc_genoprob` output by
  individual using individuals that aren't present in the data.

- Fix bug in `est_map` where it was producing `NaN`s in some cases.


## qtl2geno 0.5-17 (2017-04-17)

### Bug fixes

- `read_cross2` now unzips a `.zip` file to a separate directory, to
  avoid possibility of clashing of multiple sets of files.

- `read_cross2` will now ignore any JSON or YAML files in the `.zip`
  file that have the pattern `__MACOSX/._*`.

- `read_cross2` will stop with an error if a `.zip` file contains
  multiple JSON or multiple YAML files. If there's both a YAML and a
  JSON file, the YAML file is used and a warning is issued.

## Minor changes

- `est_map` now gives a warning if it reaches the maximum number of
  iterations without converging.


## qtl2geno 0.5-16 (2017-04-16)

### New features

- Implemented new cross types `"risib4"`, `"risib8"`, and `"magic19"`.
  The `"risib8"` cross type corresponds to the Collaborative Cross.
  The `"magic19"` cross type corresponds to the 19-way Arabidopsis
  MAGIC lines of
  [Kover et al (2009) PLOS Genetics 5:e1000551](https://doi.org/10.1371/journal.pgen.1000551).


## qtl2geno 0.5-15 (2017-04-14)

### New features

- Added argument `lowmem` to `est_map`; default is `FALSE`, which
  corresponds to a new implementation that uses more memory but is
  considerably faster.

- Added function `find_map_gaps` for identifying larger inter-marker gaps
  in a genetic map.

- Added function `calc_geno_freq` for calculating genotype
  frequencies, by individual or by marker (from the multipoint
  genotype probabilities returned by `calc_genoprob`).


## qtl2geno 0.5-14 (2017-04-05)

### New features

- Implemented new cross types `"riself4"`, `"riself8"`, and
  `"riself16"`, for multi-way MAGIC populations (multi-way RILs by
  selfing).

## Bug fixes

- Fixed problem in `read_cross2` in the case that data has a physical
  map but not a genetic map.

## Minor changes

- Added argument 'overwrite' to `write_control_file`; if `TRUE`,
  overwrite the file, if it's present. (Previously, you were always
  forced to first remove it.)


## qtl2geno 0.5-13 (2017-04-03)

### New features

- Added function `ind_ids_covar` to grab individual IDs from the
  covariate data.

- `ind_ids()` now return individuals that are in any of geno, pheno, covar.

### Bug fixes

- `subset_cross2()` now deals properly with the case that chromosome
  or individual IDs are not found in cross object, and deals with the
  case that geno and pheno (and covar) have different individuals.


## qtl2geno 0.5-12 (2017-03-30)

### New features

- Added functions `count_xo` and `locate_xo` for getting estimates of
  the number of crossovers on each chromosome in each individual, and
  of their locations.

- Added `compare_geno` for comparing raw genotypes between pairs of
  individuals (to look for possible sample duplicates).

- Added `calc_errorlod` to help identify potential genotyping errors
  (and problem markers or individuals).


## qtl2geno 0.5-9 (2017-03-23)

### Minor changes

- Made various small improvements to the handling of problems in the
  input files.

- Small changes to better handle genotype probabilities that are in
  the [qtl2feather](https://github.com/byandell/qtl2feather) format.


## qtl2geno 0.5-8 (2017-03-21)

### New features

- Added internal functions `dim.calc_genoprob` and
  `dimnames.calc_genoprob`, from
  [Brian Yandell](https://github.com/byandell), for use with
  [qtl2feather](https://github.com/byandell/qtl2feather), which uses
  [feather](https://github.com/wesm/feather) to store genotype
  probabilities in a file (to save memory).

- In precess of revising various functions to use qtl2feather,
  particularly in grabbing dimnames (with the above functions), but
  also to avoid `seq(along=genoprobs)` and instead use
  `seq_len(length(genoprobs))`.


## qtl2geno 0.5-7 (2017-03-21)

### New features

- Removed the distinction between "lines" and "individuals", and the
  `linemap` component in the input that connected them.
  (While for RILs like the Collaborative Cross, we may want to work
  with individual-level phenotypes, it seems best to deal with that
  outside of the cross object.)

- Removed the functions `n_lines()` and `line_ids()`. Added some
  functions:

    - `n_ind_geno()` for number of genotyped individuals, and
      `ind_ids_geno()` to get their IDs.
    - `n_ind_pheno()` for number of phenotyped individuals, and
      `ind_ids_pheno()` to get their IDs.
    - `n_ind_gnp()` for number of individuals with *both* genotypes
      and phenotypes, and `ind_ids_gnp()` to get their IDs.

- Also, `n_ind()` and `ind_ids()` now return the total number of
  individuals, across both genotypes and phenotypes.


## qtl2geno 0.5-6 (2017-03-17)

### New features

- Added a function `find_ibd_segments` that takes genotypes for a set
  of inbred strains and searches for segments where strain pairs look
  to be IBD.


## qtl2geno 0.5-5 (2017-03-13)

### New features

- Refactored to simplify the main data structures for genetic map and
  genotype probabilities. `calc_genoprob` now needs you to provide a
  pseudomarker map (created with `insert_pseudomarkers`).

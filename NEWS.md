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

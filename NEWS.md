## qtl2geno 0.5-11 (2017-03-29)

### New feature

- Added `compare_geno` for comparing raw genotypes between pairs of
  individuals (to look for possible sample duplicates).


## qtl2geno 0.5-10 (2017-03-27)

### New feature

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

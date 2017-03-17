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

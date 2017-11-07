## qtl2plot 0.5-17 (2017-11-07)

### New features

- Added function `plot_genoprob()` for plotting a heat map of the
  genotype probabilities for a single individual on a single
  chromosome.


## qtl2plot 0.5-16 (2017-10-17)

### New features

- Added function `plot_pxg()` for plotting phenotypes versus genotype
  at a putative QTL.


## qtl2plot 0.5-13 (2017-08-04)

### Minor changes

- Added arguments to plot_genes to make it more flexible; in
  particular, to allow differences in column names and to allow start
  and stop positions to be in either bp or Mbp.


## qtl2plot 0.5-12 (2017-07-31)

### New features

- Added `plot_peaks()`, for plotting the locations (and optionally
  confidence intervals) of inferred QTL.

- Added `plot_onegeno()` for plotting the genome-wide genotypes for
  one individual.


## qtl2plot 0.5-10 (2017-07-11)

### Bug fixes

- Ensure that `scan1()` output and the map are aligned.


## qtl2plot 0.5-7 (2017-06-05)

### Minor changes

- Revised installation instructions.


## qtl2plot 0.5-6 (2017-04-29)

### Minor changes

- Trap cases of the input `map` being `NULL`. This happens to me
  particularly when I try `some_cross$pmap` but the cross object
  doesn't contain a physical map.


## qtl2plot 0.5-4 (2017-03-11)

### Bug fixes

- Fix bug in `plot_scan1`: need to import `qtl2scan::subset_scan1`.


## qtl2plot 0.5-3 (2017-03-07)

### New features

- Refactored to deal with changes in data structures in
  [qtl2geno](https://github.com/rqtl/qtl2geno) and
  [qtl2scan](https://github.com/rqtl/qtl2scan).

- Most functions now need you to provide the map of
  markers/pseudomarkers, produced by `insert_pseudomarkers`.

- The `plot_snpasso` function requires the `snpinfo` as supplemented
  by `index_snps` (that is, with the information on which groups of
  SNPs are equivalent).

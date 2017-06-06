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

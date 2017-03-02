## qtl2plot 0.5-2

### New features

- Refactored to deal with changes in data structures in
  [qtl2geno](https://github.com/rqtl/qtl2geno) and
  [qtl2scan](https://github.com/rqtl/qtl2scan).

- Most functions now need you to provide the map of
  markers/pseudomarkers, produced by `insert_pseudomarkers`.

- The `plot_snpasso` function requires the `snpinfo` as supplemented
  by `index_snps` (that is, with the information on which groups of
  SNPs are equivalent).

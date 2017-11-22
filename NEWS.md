## qtl2plot 0.5-18 (2017-11-21)

### New features

- `plot_coef()` now takes an argument `legend` to control addition of
  a legend. If NULL, no legend is included. Otherwise it indicates the
  placement of the legend (for example `"topright"` or `"bottomleft"`.
  There are hidden arguments `legend_lab` for labels and
  `legend_ncol` for the number of columns.

- `plot_snpasso()` now takes a `minlod` argument; points with LOD less
  than this value are omitted.

- New (hidden) graphics paramater for `plot_scan1()`, `altcol`, for
  having colors alternate between chromosomes, for the case of a
  Manhattan plot particular of SNP association results.

### Minor changes

- Fixed `plot_snpasso()` so it can handle results for multiple
  chromosomes.


## qtl2plot 0.5-17 (2017-11-08)

### New features

- Added function `plot_genoprob()` for plotting a heat map of the
  genotype probabilities for a single individual on a single
  chromosome.

- Added function `plot_genoprobcomp()` for plotting a bivariate heat
  map comparing two sets of genotype probabilities for a single
  individual on a single chromosome.

- Added `swap_axes` argument to `plot_onegeno()` (to have the option of
  chromosomes running horizontally rather than vertically) and to
  `plot_pxg()` (to have the option of having genotypes on y-axis and
  phenotypes horizontally).

### Minor changes

- Throughout, changed arguments like `hlines.col` to use an underscore
  rather than a period: `hlines_col`.

- Hid a bunch of graphics parameters, to be passed via `...`, so that
  the function definitions aren't so cluttered.

- For `plot_snpasso()`, changed arguments `drop.hilit` and `col.hilit`
  to `drop_hilit` and `col_hilit` for consistency with other
  functions' argument names.


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

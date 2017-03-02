## qtl2scan 0.5-1

### New features

- Refactored to simplify the main data structures for `scan1`,
  `scan1coef`, and `scan1blup` output, and to deal with the
  refactoring of data structures in [qtl2geno](https://github.com/rqtl/qtl2geno).
  Functions like `max_scan1`, `find_peaks`, `lod_int`, and `bayes_int`
  now need you to provide a map. Similarly, to subset `scan1` results
  by chromosome, you need to provide a map to the subsetting function.

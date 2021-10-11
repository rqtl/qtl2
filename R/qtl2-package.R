#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @useDynLib qtl2, .registration=TRUE
#'
#' @section Vignettes:
#' * [user guide](https://kbroman.org/qtl2/assets/vignettes/user_guide.html)
#' * [categorized list of functions in R/qtl2](https://kbroman.org/qtl2/pages/rqtl2_functions.html)
#' * [input file formats](https://kbroman.org/qtl2/assets/vignettes/input_files.html)
#' * [using qtl2fst for on-disk genotype probabilities](https://kbroman.org/qtl2/assets/vignettes/qtl2fst.html)
#' * [preparing DO mouse data for R/qtl2](https://kbroman.org/qtl2/pages/prep_do_data.html)
#' * [genotype diagnostics for diversity outbred mice](https://kbroman.org/qtl2/assets/vignettes/do_diagnostics.html)
#' * [identifying sample mix-ups in diversity outbred mice](https://kbroman.org/qtl2/assets/vignettes/do_mixups.html)
#' * [differences between R/qtl and R/qtl2](https://kbroman.org/qtl2/assets/vignettes/rqtl_diff.html)
#' * [developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html)
#' * [Tutorial on R/qtl2](https://smcclatchy.github.io/mapping/) by
#'   [Susan McClatchy](https://github.com/smcclatchy) and [Dan Gatti](https://github.com/dmgatti)
#'
#' @section Related packages:
#' * [qtl2convert](https://github.com/rqtl/qtl2convert), for converting
#'   data among the R/qtl2, DOQTL, and R/qtl formats
#' * [qtl2fst](https://github.com/rqtl/qtl2fst), for storing genotype
#'   probabilities on disk
#' * [qtl2ggplot](https://github.com/byandell/qtl2ggplot),
#'   for [ggplot2](https://ggplot2.tidyverse.org/)-based data visualizations
#'
"_PACKAGE"

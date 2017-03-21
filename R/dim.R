#' @export
dim.calc_genoprob <- function(x) {
  sapply(x, dim)
}
#' @export
dimnames.calc_genoprob <- function(x) {
  dnames <- lapply(x, dimnames)
  list(ind = dnames[[1]][[1]],
       gen = lapply(dnames, function(x) x[[2]]),
       mar = lapply(dnames, function(x) x[[3]]))
}
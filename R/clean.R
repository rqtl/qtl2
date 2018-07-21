#' Clean an object
#'
#' Clean an object by removing messy values
#'
#' @md
#'
#' @param object Object to be cleaned
#' @param ... Other arguments
#' @return Input object with messy values cleaned up
#'
#' @export
#' @seealso [clean.scan1()], [clean.calc_genoprob()]
clean <-
    function(object, ...)
{
    UseMethod("clean")
}

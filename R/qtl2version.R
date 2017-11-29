#' Installed version of R/qtl2
#'
#' Get installed version of R/qtl2
#'
#' @return A character string with the installed version of the R/qtl2 package.
#'
#' @examples
#' qtl2version()
#'
#' @export
qtl2version <-
    function()
{
    version <- unlist(packageVersion("qtl2"))

    # make it like #.#-#
    paste(c(version, ".", "-")[c(1,4,2,5,3)], collapse="")
}

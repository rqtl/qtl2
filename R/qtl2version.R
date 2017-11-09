#' List R/qtl2 package versions.
#'
#' List versions of all installed R/qtl2 packages, as shown in the
#' results of [utils::installed.packages()].
#'
#' @md
#'
#' @return Data frame with two columns: package name and version.
#'
#' @examples
#' qtl2version()
#'
#' @importFrom utils installed.packages
#' @export

qtl2version <-
    function()
{
    # grab all installed packages
    z <- installed.packages()

    # pull out those whose names start with "qtl2"
    qtl2 <- grep("^qtl2", z[,"Package"])

    # pull out the key columns, as a data frame
    result <- as.data.frame(z[qtl2,c("Package", "Version"),drop=FALSE],
                            stringsAsFactors=FALSE)

    # reorder to put core packages first
    core <- c("qtl2", "qtl2convert", "qtl2db", "qtl2geno", "qtl2plot", "qtl2scan", "qtl2bioc")
    result <- result[order(result$Package %in% core, decreasing=TRUE),,drop=FALSE]

    # make the rownames numbers
    rownames(result) <- 1:nrow(result)

    result
}

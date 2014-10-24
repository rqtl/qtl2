# probs_to_grid
#' Subset genotype probability array to pseudomarkers on a grid
#'
#' Subset genotype probability array (from \code{\link{calc_genoprob}}
#' to a grid of pseudomarkers along each chromosome.
#'
#' @param probs List of 3d arrays, as output from
#' \code{\link{calc_genoprob}} with \code{stepwidth="fixed"}.
#'
#' @return Same list as input, but subset to just include
#' pseudomarkers along a grid.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' sapply(probs, dim)
#' probs_sub <- probs_to_grid(probs)
#' sapply(probs_sub, dim)

probs_to_grid <-
    function(probs)
{
    if(!("map" %in% names(attributes(probs))))
        stop("probs has no map attribute")

    # check stepwidth arg
    map <- attr(probs, "map")
    stepwidth <- vapply(map, attr, "", "stepwidth")
    if(!all(stepwidth=="fixed"))
        stop('probs needs to be result of calc_genoprob with stepwidth="fixed"')

    for(i in seq(along=probs)) {
        # grab grid vector
        grid <- attr(map[[i]], "grid")
        if(is.null(grid))
            stop("grid attribute not found for chr ", names(probs)[i])

        # subset probs
        probs[[i]] <- probs[[i]][,grid,]
    }

    # attribute indicating that the thing has been subsetted
    attr(probs, "subset") <- TRUE

    probs
}

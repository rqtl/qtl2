# probs_to_grid
#' Subset genotype probability array to pseudomarkers on a grid
#'
#' Subset genotype probability array (from \code{\link{calc_genoprob}}
#' to a grid of pseudomarkers along each chromosome.
#'
#' @param probs Genotype probabilities as output from
#' \code{\link{calc_genoprob}} with \code{stepwidth="fixed"}.
#'
#' @return Same list as input, but subset to just include
#' pseudomarkers along a grid. The map attribute is similarly subset.
#'
#' @details This only works if \code{\link{calc_genoprob}} was run
#' with \code{stepwidth="fixed"}, so that the genotype probabilities
#' were calculated at a grid of markers/pseudomarkers. When this is
#' the case, we omit all but the probabilities on this grid.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' sapply(probs$probs, dim)
#' probs_sub <- probs_to_grid(probs)
#' sapply(probs_sub$probs, dim)

probs_to_grid <-
    function(probs)
{
    map <- probs$map
    if(is.null(probs$map))
        stop("probs has no map attribute")
    # check stepwidth arg
    stepwidth <- probs$stepwidth
    if(is.null(stepwidth) || stepwidth != "fixed")
        stop('probs needs to be result of calc_genoprob with stepwidth="fixed"')

    grid <- probs$grid
    if(is.null(probs$grid))
        stop("probs has no grid attribute")

    chrID <- names(probs$probs)
    for(i in seq(along=chrID)) {
        # grab grid vector
        if(is.null(grid[[i]])) {
            stop("grid not found for chr ", chrID[i])
        }
        if(length(grid[[i]]) != length(map[[i]])) {
            stop("length(grid) [", length(grid[[i]]), "] != length(map) [",
                 length(map[[i]]), "] for chr ", chrID[i])
        }

        # subset probs
        if(!all(grid[[i]])) {
            if(length(grid[[i]]) != dim(probs$probs[[i]])[3])
                stop("length(grid) [", length(grid[[i]]), "] != ncol(probs) [",
                     ncol(probs$probs[[i]]), "] for chr ", chrID[i])
            probs$probs[[i]] <- probs$probs[[i]][,,grid[[i]],drop=FALSE]
        }
    }

    probs$map <- map_to_grid(map, grid)
    probs$grid <- NULL # don't need this anymore

    probs
}

# subset a map object to grid
#
# input is a list; attributes include "grid"
map_to_grid <-
    function(map, grid)
{
    if(is.null(grid)) stop("grid is NULL")

    for(i in seq(along=map)) {
        if(is.null(grid[[i]]) || all(grid[[i]])) next

        # subset map
        map[[i]] <- map[[i]][grid[[i]]]
    }

    map
}

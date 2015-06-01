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
        if(!all(grid)) {
            if(length(grid) != ncol(probs[[i]]))
                stop("length(grid) (", length(grid), ") != ncol(probs) (",
                     ncol(probs[[i]]), ") for chr ", names(probs)[i])
            probs[[i]] <- probs[[i]][,grid,]
        }
    }

    attr(probs, "map") <- map_to_grid(map)

    probs
}

# subset a map object to grid
#
# input is a list; attributes include "grid"
map_to_grid <-
    function(map)
{
    for(i in seq(along=map)) {
        mapat <- attributes(map[[i]])
        grid <- mapat$grid
        if(is.null(grid) || all(grid)) next

        # subset map
        map[[i]] <- map[[i]][grid]

        mapat_ignore <- "names"
        mapat_subset <- c("index", "grid")
        for(att in names(mapat)) {
            if(att %in% mapat_ignore) next
            if(att %in% mapat_subset)
                attr(map[[i]], att) <- mapat[[att]][grid]
            else
                attr(map[[i]], att) <- mapat[[att]]
        }
    }

    map
}


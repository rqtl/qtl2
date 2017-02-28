# probs_to_grid
#' Subset genotype probability array to pseudomarkers on a grid
#'
#' Subset genotype probability array (from \code{\link{calc_genoprob}}
#' to a grid of pseudomarkers along each chromosome.
#'
#' @param probs Genotype probabilities as output from
#' \code{\link{calc_genoprob}} with \code{stepwidth="fixed"}.
#'
#' @param grid List of logical vectors that indicate which positions
#' are on the grid and should be retained.
#'
#' @return Same list as input, but subset to just include
#' pseudomarkers along a grid. The map attribute is similarly subset.
#'
#' @details This only works if \code{\link{calc_genoprob}} was run
#' with \code{stepwidth="fixed"}, so that the genotype
#' probabilities were calculated at a grid of
#' markers/pseudomarkers. When this is the case, we omit all but
#' the probabilities on this grid. Use \code{\link{calc_grid}} to
#' find the grid positions.
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{calc_grid}}, \code{\link{map_to_grid}}
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' map_w_pmar <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map_w_pmar, error_prob=0.002)
#' sapply(probs, dim)
#' grid <- calc_grid(grav2$gmap, step=1)
#' probs_sub <- probs_to_grid(probs, grid)
#' sapply(probs_sub, dim)

probs_to_grid <-
    function(probs, grid)
{
    chrID <- names(probs)
    for(i in seq(along=chrID)) {
        # grab grid vector
        if(is.null(grid[[i]]) || all(grid[[i]])) next

        # subset probs
        if(length(grid[[i]]) != dim(probs[[i]])[3])
            stop("length(grid) [", length(grid[[i]]), "] != ncol(probs) [",
                 ncol(probs[[i]]), "] for chr ", chrID[i])
        probs[[i]] <- probs[[i]][,,grid[[i]],drop=FALSE]
    }

    probs
}

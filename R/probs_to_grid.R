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
#' map_w_pmar <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map_w_pmar, error_prob=0.002)
#' sapply(probs, dim)
#' probs_sub <- probs_to_grid(probs)
#' sapply(probs_sub, dim)

probs_to_grid <-
    function(probs)
{
    grid <- attr(probs, "grid")
    if(is.null(grid)) {
        warning("probs has no grid attribute")
        return(probs)
    }

    chrID <- names(probs)
    for(i in seq(along=chrID)) {
        # grab grid vector
        if(is.null(grid[[i]])) {
            stop("grid not found for chr ", chrID[i])
        }
        # subset probs
        if(!all(grid[[i]])) {
            if(length(grid[[i]]) != dim(probs[[i]])[3])
                stop("length(grid) [", length(grid[[i]]), "] != ncol(probs) [",
                     ncol(probs[[i]]), "] for chr ", chrID[i])
            probs[[i]] <- probs[[i]][,,grid[[i]],drop=FALSE]
        }
    }

    attr(probs, "grid") <- NULL # don't need this anymore

    probs
}

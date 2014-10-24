# calc_genetic_sim
#' Calculate genetic similarity among individuals
#'
#' Calculate genetic similarity among individuals from conditional genotype probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as calculated from \code{\link{calc_genoprob}}.
#' @param use_grid_only If \code{TRUE} and \code{probs} were calculated with
#' \code{stepwidth="fixed"}, reduce them to the grid using
#' \code{\link{probs_to_grid}}.
#' @param omit_x If \code{TRUE}, only use the autosomes.
#' @param quiet IF \code{FALSE}, print messages on progress.
#'
#' @return A matrix of proportion of matching genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' sim <- calc_genetic_sim(probs)

calc_genetic_sim <-
    function(probs, use_grid_only=TRUE, omit_x=TRUE, quiet=TRUE)
{
    n_ind <- nrow(probs[[1]])
    ind_names <- rownames(probs[[1]])
    result <- matrix(0, ncol=n_ind, nrow=n_ind)
    dimnames(result) <- list(ind_names, ind_names)

    if(omit_x) chr <- which(!attr(probs, "is_x_chr"))
    else chr <- seq(along=probs)

    subsetted <- attr(probs, "subset")
    stepwidth <- attr(attr(probs, "map")[[1]], "stepwidth")
    if(use_grid_only &&
       (is.null(subsetted) || !subsetted) &&
       !is.null(stepwidth) && stepwidth=="fixed") {
        if(!quiet) message(" - Reducing probabilities to grid")
        probs <- probs_to_grid(probs)
    }

    ## FIX ME: should reduce to allele counts

    tot_pos <- 0
    for(i in chr) {
        if(!quiet) message(" - Chr ", names(probs)[i])
        result <- result + .calc_genetic_sim(aperm(probs[[i]], c(2,3,1)))
        tot_pos <- tot_pos + ncol(probs[[i]])
    }

    result/tot_pos
}

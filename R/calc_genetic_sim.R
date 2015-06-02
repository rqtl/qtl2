# calc_genetic_sim
#' Calculate genetic similarity among individuals
#'
#' Calculate genetic similarity among individuals from conditional genotype probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as
#' calculated from \code{\link{calc_genoprob}}.
#' @param use_grid_only If \code{TRUE} and \code{probs} were calculated with
#' \code{stepwidth="fixed"}, reduce them to the grid using
#' \code{\link{probs_to_grid}}.
#' @param omit_x If \code{TRUE}, only use the autosomes.
#' @param use_allele_probs If \code{TRUE}, assess similarity with
#' allele probabilities (that is, first run
#' \code{\link{genoprob_to_alleleprob}}); otherwise use the genotype
#' probabilities.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param n_cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#'
#' @return A matrix of proportion of matching alleles.
#'
#' @details If \code{use_allele_probs=TRUE} (the default), we first
#' convert the genotype probabilities are converted to allele
#' probabilities (using \code{\link{genoprob_to_alleleprob}}). This is
#' recommended, as then the result is like a empirical kinship
#' coefficient (e.g., the expected value for an intercross is 1/2;
#' using genotype probabilities, the expected value is 3/8).
#'
#' We then calculate
#' \eqn{\sum_{kl}(p_{ikl} p_{jkl})}{sum_kl (p_ikl p_jkl)}
#' where \eqn{k} = position, \eqn{l} = allele, and \eqn{i,j}
#' are two individuals.
#'
#' For crosses with just two possible genotypes (e.g., backcross), we
#' don't convert to allele probabilities but just use the original
#' genotype probabilities.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' sim <- calc_genetic_sim(probs)

calc_genetic_sim <-
    function(probs, use_grid_only=TRUE, omit_x=TRUE,
             use_allele_probs=TRUE, quiet=TRUE, n_cores=1)
{
    n_ind <- nrow(probs[[1]])
    ind_names <- rownames(probs[[1]])
    result <- matrix(0, ncol=n_ind, nrow=n_ind)
    dimnames(result) <- list(ind_names, ind_names)

    if(omit_x) chrs <- which(!attr(probs, "is_x_chr"))
    else chrs <- seq(along=probs)

    if(n_cores==0) n_cores <- parallel::detectCores() # if 0, detect cores
    if(n_cores > 1) {
        if(!quiet) message(" - Using ", n_cores, " cores.")
        quiet <- TRUE # no more messages
    }

    stepwidth <- attr(attr(probs, "map")[[1]], "stepwidth")
    if(use_grid_only &&
       !is.null(stepwidth) && stepwidth=="fixed") {
        if(!quiet) message(" - Reducing probabilities to grid")
        probs <- probs_to_grid(probs)
    }

    # convert from genotype probabilities to allele probabilities
    if(use_allele_probs) {
        if(!quiet) message(" - converting to allele probs")
        probs <- genoprob_to_alleleprob(probs, quiet=quiet, n_cores=n_cores)
    }

    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        .calc_genetic_sim(aperm(probs[[chr]], c(3,2,1))) # convert to pos x gen x ind
    }

    if(n_cores<=1) { # no parallel processing
        for(chr in chrs)
            result <- result + by_chr_func(chr)
    }
    else {
        if(Sys.info()[1] == "Windows") { # Windows doesn't suport mclapply
            cl <- parallel::makeCluster(n_cores)
            on.exit(parallel::stopCluster(cl))
            by_chr_res <- parallel::clusterApply(cl, chrs, by_chr_func)
        }
        else {
            by_chr_res <- parallel::mclapply(chrs, by_chr_func, mc.cores=n_cores)
        }
        for(chr in seq(along=by_chr_res))
            result <- result + by_chr_res[[chr]]
    }

    tot_pos <- sum(vapply(probs, function(a) dim(a)[3], 0)[chrs])
    result/tot_pos
}

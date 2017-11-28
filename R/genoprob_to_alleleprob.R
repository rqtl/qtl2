# genoprob_to_alleleprob
#' Convert genotype probabilities to allele probabilities
#'
#' Reduce genotype probabilities (as calculated by
#' [calc_genoprob()]) to allele probabilities.
#'
#' @md
#'
#' @param probs Genotype probabilities, as calculated from
#' [calc_genoprob()].
#' @param quiet IF `FALSE`, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return The `probs` input with probabilities
#' collapsed to alleles rather than genotypes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' gmap_w_pmar <- insert_pseudomarkers(iron, step=1)
#' probs <- calc_genoprob(iron, gmap_w_pmar, error_prob=0.002)
#' allele_probs <- genoprob_to_alleleprob(probs)

genoprob_to_alleleprob <-
    function(probs, quiet=TRUE, cores=1)
{
    # already converted?
    ap <- attr(probs, "alleleprobs")
    if(!is.null(ap) && ap) return(probs)

    is_x_chr <- attr(probs, "is_x_chr")

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # no more messages
    }

    probs_attr <- attributes(probs)

    # alleles attribute?
    alleles <- probs_attr$alleles
    if(is.null(alleles))
        warning("probs has no alleles attribute; guessing allele codes.")

    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        result <- aperm(.genoprob_to_alleleprob(attr(probs, "crosstype"),
                                                aperm(probs[[chr]], c(2, 1, 3)), # reorg -> geno x ind x pos
                                                is_x_chr[chr]),
                        c(2, 1, 3)) # reorg back to ind x geno x pos

        # allele names
        dn <- dimnames(probs)
        if(is.null(alleles) || length(alleles) < ncol(result)) {
            alleles <- assign_allele_codes(ncol(result), dn[[2]][[chr]])
        }
        dn[[2]][[chr]] <- alleles
        dimnames(result) <- list(dn[[1]], dn[[2]][[chr]], dn[[3]][[chr]])

        result
    }

    chrID <- names(probs)
    chrs <- seq(along=chrID)

    probs <- cluster_lapply(cores, chrs, by_chr_func) # if cores==1, this uses lapply()
    names(probs) <- chrID

    attr(probs, "crosstype") <- probs_attr$crosstype
    attr(probs, "is_x_chr") <- probs_attr$is_x_chr
    attr(probs, "alleles") <- probs_attr$alleles
    attr(probs, "alleleprobs") <- TRUE
    class(probs) <- c("calc_genoprob", "list")

    probs
}

# calc_kinship
#' Calculate kinship matrix
#'
#' Calculate genetic similarity among individuals (kinship matrix)
#' from conditional genotype probabilities.
#'
#' @param probs List of three-dimensional arrays of probabilities, as
#' calculated from \code{\link{calc_genoprob}}.
#' @param type Indicates whether to calculate the overall kinship
#' (\code{"overall"}, using all chromosomes), the kinship matrix
#' leaving out one chromosome at a time (\code{"loco"}), or the
#' kinship matrix for each chromosome (\code{"chr"}).
#' @param use_grid_only If \code{TRUE} and \code{probs} were calculated with
#' \code{stepwidth="fixed"}, reduce them to the grid using
#' \code{\link{probs_to_grid}}.
#' @param omit_x If \code{TRUE}, only use the autosomes; ignored when
#' \code{type="chr"}.
#' @param use_allele_probs If \code{TRUE}, assess similarity with
#' allele probabilities (that is, first run
#' \code{\link{genoprob_to_alleleprob}}); otherwise use the genotype
#' probabilities.
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return If \code{type="overall"} (the default), a matrix of
#' proportion of matching alleles. Otherwise a list with one matrix
#' per chromosome.
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
#' sim <- calc_kinship(probs)

calc_kinship <-
    function(probs, type=c("overall", "loco", "chr"),
             use_grid_only=TRUE, omit_x=TRUE,
             use_allele_probs=TRUE, quiet=TRUE, cores=1)
{
    type <- match.arg(type)

    allchr <- names(probs)
    if(omit_x && type != "chr") chrs <- which(!attr(probs, "is_x_chr"))
    else chrs <- seq(along=allchr)

    stepwidth <- attr(attr(probs, "map")[[1]], "stepwidth")
    if(use_grid_only &&
       !is.null(stepwidth) && stepwidth=="fixed") {
        if(!quiet) message(" - Reducing probabilities to grid")
        probs <- probs_to_grid(probs)
    }

    # convert from genotype probabilities to allele probabilities
    if(use_allele_probs) {
        if(!quiet) message(" - converting to allele probs")
        probs <- genoprob_to_alleleprob(probs, quiet=quiet, cores=cores)
    }

    if(type=="overall")
        return(calc_kinship_overall(probs, chrs=chrs, quiet=quiet, cores=cores))
    else if(type=="chr")
        return(calc_kinship_bychr(probs, chrs=chrs, scale=TRUE, quiet=quiet, cores=cores))

    # otherwise LOCO (leave one chromosome out)
    result <- calc_kinship_bychr(probs, chrs=chrs, scale=FALSE, quiet=quiet, cores=cores)
    kinship_bychr2loco(result, allchr)
}

# calculate an overall kinship matrix
calc_kinship_overall <-
    function(probs, chrs, quiet=TRUE, cores=1)
{
    n_ind <- nrow(probs[[1]])
    ind_names <- rownames(probs[[1]])

    result <- matrix(0, nrow=n_ind, ncol=n_ind)
    dimnames(result) <- list(ind_names, ind_names)

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # function that does the work
    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])
        pr <- aperm(probs[[chr]], c(3,2,1)) # convert to pos x gen x ind
        .calc_kinship(pr)
    }

    # run and combine results
    if(n_cores(cores) == 1) {
        for(chr in chrs)
            result <- result + by_chr_func(chr)
    }
    else {
        by_chr_res <- cluster_lapply(cores, chrs, by_chr_func)
        for(i in seq(along=by_chr_res))
            result <- result + by_chr_res[[i]]
    }

    tot_pos <- sum(vapply(probs, function(a) dim(a)[3], 0)[chrs])
    result <- result/tot_pos
    attr(result, "n_pos") <- tot_pos
    result
}

# calculate kinship for each chromosome
calc_kinship_bychr <-
    function(probs, chrs, scale=TRUE, quiet=TRUE, cores=1)
{
    n_ind <- nrow(probs[[1]])
    ind_names <- rownames(probs[[1]])

    # set up cluster and set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        cat(" - Using", n_cores(cores), "cores\n")
        quiet <- TRUE # make the rest quiet
    }

    # function that does the work
    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", names(probs)[chr])

        n_pos <- dim(probs[[chr]])[3]

        # aperm converts to pos x gen x ind
        pr <- aperm(probs[[chr]], c(3,2,1))
        result <- .calc_kinship(pr)
        if(scale) result <- result/n_pos

        attr(result, "n_pos") <- n_pos
        dimnames(result) <- list(ind_names, ind_names)
        result
    }

    # run and combine results
    result <- cluster_lapply(cores, chrs, by_chr_func)

    names(result) <- names(probs)[chrs]
    result
}

# use kinship for each chromosome
# to calculate kinship leaving one chromosome out at a time
kinship_bychr2loco <-
    function(kinship, allchr)
{
    # sum over chromosomes
    overall <- kinship[[1]]
    tot_pos <- attr(kinship[[1]], "n_pos")
    if(length(kinship) > 1) {
        for(i in 2:length(kinship)) {
            overall <- overall + kinship[[i]]
            tot_pos <- tot_pos + attr(kinship[[i]], "n_pos")
        }
    }
    attr(overall, "n_pos") <- tot_pos

    for(chr in allchr) {
        if(chr %in% names(kinship)) {
            n_pos <- attr(kinship[[chr]], "n_pos")
            kinship[[chr]] <- (overall - kinship[[chr]])/(tot_pos - n_pos)
            attr(kinship[[chr]], "n_pos") <- tot_pos - n_pos
        } else {
            kinship <- c(kinship, list(overall/tot_pos))
            names(kinship)[length(kinship)] <- chr
        }
    }

    # make sure it's in the right order
    kinship[allchr]
}

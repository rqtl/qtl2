# calc_errorlod
#' Calculate genotyping error LOD scores
#'
#' Use the genotype probabilities calculated with
#' [calc_genoprob()] to calculate genotyping error LOD
#' scores, to help identify potential genotyping errors (and problem
#' markers and/or individuals).
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param probs Genotype probabilities as calculated from [calc_genoprob()].
#' @param quiet If `FALSE`, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list of matrices of genotyping error LOD scores. Each
#'     matrix corresponds to a chromosome and is arranged as
#'     individuals x markers.
#'
#' @details
#'   Let \eqn{O_k}{O[k]} denote the observed marker genotype at position
#'  \eqn{k}, and \eqn{g_k}{g[k]} denote the corresponding true underlying
#'  genotype.
#'
#'  Following Lincoln and Lander (1992), we calculate
#'  LOD = \eqn{log_{10} [ Pr(O_k | g_k = O_k) / Pr(O_k | g_k \ne O_K) ]}{
#'  log10[ Pr(O[k] | g[k] = O[k]) / Pr(O[k] | g[k] != O[k]) ]}
#'
#' @references
#' Lincoln SE, Lander ES (1992) Systematic detection of errors in genetic linkage data. Genomics 14:604--610.
#'
#' @export
#' @keywords utilities
#' @seealso [calc_genoprob()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' probs <- calc_genoprob(iron, error_prob=0.002, map_function="c-f")
#' errorlod <- calc_errorlod(iron, probs)
#'
#' # combine into one matrix
#' errorlod <- do.call("cbind", errorlod)

calc_errorlod <-
function(cross, probs, quiet=TRUE, cores=1)
{
    # check inputs
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    alleleprobs <- attr(probs, "alleleprobs")
    if(!is.null(alleleprobs) && alleleprobs)
        stop("probs must contain full genotype probabilities, not haplotype/allele probabilities")

    if(length(cross$geno) != length(probs) ||
       !all(names(cross$geno) == names(probs))) {
        chr1 <- names(cross$geno)
        chr2 <- names(probs)
        chr <- find_common_ids(chr1, chr2)
        if(length(chr)==0)
            stop("cross and probs have no chromosomes in common")
        cross <- cross[,chr]
        probs <- subset(probs, chr=chr)
    }

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # deal with missing information
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    dn <- dimnames(probs)

    # align individuals
    if(length(ind) != length(dn[[1]]) ||
       !all(ind == dn[[1]])) {
        ind <- find_common_ids(ind, dn[[1]])
        if(length(ind)==0)
            stop("No individuals in common between cross and probs")
        cross <- cross[ind,]
    }

    by_group_func <- function(i) {
        mn <- colnames(cross$geno[[chr]])
        if(!all(mn %in% dn[[3]][[chr]]))
            stop("Some markers in cross are not in probs on chr ", chrnames[chr])
        pr <- aperm(probs[[chr]][ind,,mn,drop=FALSE], c(2,1,3)) # genotype, ind, marker
        errorlod <- t(.calc_errorlod(cross$crosstype, pr[,group[[i]],,drop=FALSE],
                                     t(cross$geno[[chr]][group[[i]],,drop=FALSE]),
                                     founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[group[[i]][1]],
                                     t(cross$cross_info[group[[i]][1],])))

    }

    # split individuals into groups with common sex and cross_info
    sex_crossinfo <- paste(cross$is_female, apply(cross$cross_info, 1, paste, collapse=":"), sep=":")
    group <- split(seq(along=sex_crossinfo), sex_crossinfo)
    names(group) <- NULL
    nc <- n_cores(cores)
    while(nc > length(group) && max(sapply(group, length)) > 1) { # successively split biggest group in half until there are as many groups as cores
        mx <- which.max(sapply(group, length))
        g <- group[[mx]]
        group <- c(group, list(g[seq(1, length(g), by=2)]))
        group[[mx]] <- g[seq(2, length(g), by=2)]
    }
    groupindex <- seq(along=group)

    errorlod <- vector("list", length(probs))
    names(errorlod) <- chrnames

    for(chr in seq_len(length(chrnames))) {
        if(!quiet) message("Chr ", chrnames[chr])

        temp <- cluster_lapply(cores, groupindex, by_group_func)

        # paste them back together
        d <- vapply(temp, dim, c(0,0))
        nr <- sum(d[1,])
        errorlod[[chr]] <- matrix(nrow=sum(vapply(temp, nrow, 1)),
                                  ncol=ncol(temp[[1]]))
        for(i in groupindex)
            errorlod[[chr]][group[[i]],] <- temp[[i]]

        dimnames(errorlod[[chr]]) <- dimnames(cross$geno[[chr]])
    }

    errorlod
}

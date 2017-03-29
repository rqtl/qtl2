# compare_geno
#' Compare individuals' genotype data
#'
#' Count the number of matching genotypes between all pairs of
#' individuals, to look for unusually closely related individuals.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param omit_x If TRUE, only use autosomal genotypes
#' @param quiet IF \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A square matrix; diagonal is number of observed genotypes
#' for each individual; upper triangle is the number of matching
#' genotypes for each pair; lower triangle is the number of
#' markers where both of a pair were genotyped. The object is
#' given class \code{"compare_geno"}.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' cg <- compare_geno(grav2)
#' summary(cg)

compare_geno <-
    function(cross, omit_x=FALSE, quiet=TRUE, cores=1)
{
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')

    # which chr is X?
    is_x_chr <- cross$is_x_chr
    if(is.null(is_x_chr))
        is_x_chr <- rep(FALSE, length(cross$geno))

    # chr names and numeric indices
    chrnames <- names(cross$geno)
    chrnum <- seq_along(chrnames)
    if(omit_x) { # don't use the X chromosome
        chrnames <- chrnames[!is_x_chr]
        chrnum <- chrnum[!is_x_chr]
    }

    ind_names <- rownames(cross$geno[[1]])
    n_ind <- length(ind_names)

    result <- matrix(0, nrow=n_ind, ncol=n_ind)
    dimnames(result) <- list(ind_names, ind_names)

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # function that does the work
    by_chr_func <- function(chr) {
        if(!quiet) message(" - Chr ", chrnames[chr])
        .compare_geno(t(cross$geno[[chr]]))
    }

    # run and combine results
    if(n_cores(cores) == 1) {
        for(chr in chrnum)
            result <- result + by_chr_func(chr)
    }
    else {
        result_list <- cluster_lapply(cores, chrnum, by_chr_func)
        for(i in seq_along(result_list))
            result <- result + result_list[[i]]
    }

    class(result) <- c("compare_geno", "matrix")
    result
}


#' Basic summary of compare_geno object
#'
#' From results of \code{\link{compare_geno}}, show individuals with similar genotypes.
#'
#' @param object A square matrix with genotype comparisons for pairs
#'     of individuals, as output by \code{\link{compare_geno}}.
#' @param threshold Minimum proportion matches for a pair of individuals to be shown.
#' @param ... Ignored
#'
#' @return Data frame with individual pair, proportion matches, number
#'     of mismatches, number of matches, and total markers genotyped.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' cg <- compare_geno(grav2)
#' summary(cg)

summary_compare_geno <-
    function(object, threshold=0.9, ...)
{
    # proportion matching
    p <- object
    p[upper.tri(p)] <- p[upper.tri(p)] / t(p)[upper.tri(p)]
    p[lower.tri(p)] <- t(p)[lower.tri(p)]
    diag(p) <- NA

    if(sum(!is.na(p) & p >= threshold) == 0) { # no results
        result <- data.frame(ind1=character(0),
                             ind2=character(0),
                             prop_match=numeric(0),
                             n_mismatch=numeric(0),
                             n_typed=numeric(0),
                             n_match=numeric(0),
                             stringsAsFactors=FALSE)

        class(result) <- c("summary.compare_geno", "data.frame")
        attr(result, "threshold") <- threshold
        return(result)
    }

    wh <- which(!is.na(p) & p >= threshold, arr.ind=TRUE)
    wh <- wh[wh[,1] < wh[,2], , drop=FALSE]

    ind <- rownames(object)

    r <- 1:nrow(wh)
    p_match <- vapply(r, function(i) p[wh[i,1],wh[i,2]], 1)
    n_match <- vapply(r, function(i) object[wh[i,1],wh[i,2]], 1)
    n_typed <- vapply(r, function(i) object[wh[i,2],wh[i,1]], 1)

    result <- data.frame(ind1=ind[wh[,1]],
                         ind2=ind[wh[,2]],
                         prop_match=p_match,
                         n_mismatch=n_typed-n_match,
                         n_typed=n_typed,
                         n_match=n_match,
                         stringsAsFactors=FALSE)

    result <- result[order(result$prop_match, decreasing=TRUE),,drop=FALSE]

    attr(result, "threshold") <- threshold
    class(result) <- c("summary.compare_geno", "data.frame")
    result
}

#' @rdname summary_compare_geno
#' @export
summary.compare_geno <- summary_compare_geno

#' @rdname summary_compare_geno
#' @param x Results of \code{\link{summary.compare_geno}}
#' @param digits Number of digits to print
#' @export
print.summary.compare_geno <-
    function(x, digits=2, ...)
{
    if(nrow(x) == 0)
        cat("No pairs with proportion matches above ", attr(x, "threshold"), ".\n", sep="")
    else {
        print.data.frame(x, digits=digits, row.names=FALSE)
    }
    invisible(x)
}

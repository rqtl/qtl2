#' Calculate heterozygosities
#'
#' Calculate heterozygosites, by individual or by marker
#'
#' @md
#'
#' @param probs List of arrays of genotype probabilities, as
#'     calculated by [calc_genoprob()].
#' @param by Whether to summarize by individual or marker
#'
#' @param omit_x If TRUE, omit the X chromosome.
#'
#' @return
#' The result is a vector of estimated heterozygosities
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' p <- calc_genoprob(iron, err=0.002)
#'
#' # heterozygosities by individual
#' het_ind <- calc_het(p)
#'
#' # heterozygosities by marker
#' het_mar <- calc_het(p, "marker")
#'
#' @export
calc_het <-
    function(probs, by=c("individual", "marker"), omit_x=TRUE)
{
    by <- match.arg(by)
    if(is.null(probs)) stop("probs is NULL")
    if("cross2" %in% class(probs))
        stop('Input probs is a "cross2" object but should be genotype probabilities, as from calc_genoprob')

    is_x_chr <- attr(probs, "is_x_chr")
    if(any(is_x_chr) && !all(is_x_chr) && omit_x) {
        probs <- probs[,!is_x_chr]
        is_x_chr <- is_x_chr[!is_x_chr]
        omit_x <- FALSE
    }

    n_chr <- length(probs)

    # determine which columns are het
    het_col <- vector("list", n_chr)
    for(chr in seq_len(n_chr)) {
        geno <- colnames(probs[[chr]])
        a1 <- substr(geno, 1, 1)
        a2 <- substr(geno, 2, 2)
        if(is_x_chr[chr]) het_col[[chr]] <- (a1 != a2 & a2 != "Y")
        else het_col[[chr]] <- (a1 != a2)
    }

    if(by=="individual") {
        # total markers
        total_mar <- sum( vapply(probs, function(a) dim(a)[3], 1) )

        # summarize each chromosome
        result <- lapply(seq_len(n_chr), function(chr) apply(probs[[chr]][,het_col[[chr]],,drop=FALSE], 1, sum))

        if(length(result)>1) {
            for(i in seq_along(result)[-1])
                result[[1]] <- result[[1]] + result[[i]]
        }
        return(result[[1]]/total_mar)
    }

    # else: by marker
    unlist(lapply(seq_len(n_chr), function(chr) apply(probs[[chr]][,het_col[[chr]],,drop=FALSE], 3, sum)/
                                                          nrow(probs[[chr]])))
}

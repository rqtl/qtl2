#' Calculate genotype frequencies
#'
#' Calculate genotype frequencies, by individual or by marker
#'
#' @param probs List of arrays of genotype probabilities, as
#'     calculated by [calc_genoprob()].
#' @param by Whether to summarize by individual or marker
#' @param omit_x If TRUE, results are just for the autosomes. If
#'     FALSE, results are a list of length two, containing the results
#'     for the autosomes and those for the X chromosome.
#'
#' @return
#' If `omit_x=TRUE`, the result is a matrix of genotype
#' frequencies; columns are genotypes and rows are either individuals
#' or markers.
#'
#' If necessary (that is, if `omit_x=FALSE`, the data include the
#' X chromosome, and the set of genotypes on the X chromosome are
#' different than on the autosomes), the result is a list with two
#' components (for the autosomes and for the X chromosome), each being
#' a matrix of genotype frequencies.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' p <- calc_genoprob(iron, err=0.002)
#'
#' # genotype frequencies by marker
#' tab_g <- calc_geno_freq(p, "marker")
#'
#' # allele frequencies by marker
#' ap <- genoprob_to_alleleprob(p)
#' tab_a <- calc_geno_freq(ap, "marker")
#'
#' @export
calc_geno_freq <-
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

    if(any(is_x_chr) && any(!is_x_chr)) { # some autosome, some X chr
        ng <- dim(probs)[2,]
        g <- dimnames(probs)[[2]]
        if(length(unique(ng)) > 1 ||  # not all the same number of genotypes
           !all(vapply(g[-1], function(a) all(a==g[[1]]),TRUE))) { # not all the same genotypes
            return( list(A=calc_geno_freq(probs[,!is_x_chr], by, FALSE),
                         X=calc_geno_freq(probs[,is_x_chr], by, FALSE)) )
        }
    }

    # for rest, can assume that they're all one group

    if(by=="individual") {
        # total markers
        total_mar <- sum( vapply(probs, function(a) dim(a)[3], 1) )

        # summarize each chromosome
        result <- lapply(probs, apply, 1:2, sum)

        if(length(result)>1) {
            for(i in seq_along(result)[-1])
                result[[1]] <- result[[1]] + result[[i]]
        }
        return(result[[1]]/total_mar)
    }

    # else: by marker
    t(do.call("cbind", lapply(probs, apply, 2:3, mean)))
}

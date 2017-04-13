#' Calculate genotype frequencies
#'
#' Calculate genotype frequencies, by individual or by marker
#'
#' @param probs List of arrays of genotype probabilities, as
#'     calculated by \code{\link{calc_genoprob}}.
#' @param by Whether to summarize by individual or marker
#' @param omit_x If TRUE, results are just for the autosomes. If
#'     FALSE, results are a list of length two, containing the results
#'     for the autosomes and those for the X chromosome.
#'
#' @return
#' If \code{omit_x=TRUE}, the result is a matrix of genotype
#' frequencies; columns are genotypes and rows are either individuals
#' or markers.
#'
#' If \code{omit_x=FALSE} and the data include the X chromosome, the
#' result is a list with two components (for the autosomes and for the
#' X chromosome), each being a matrix of genotype frequencies.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' p <- calc_genoprob(iron, err=0.002)
#' tab <- calc_geno_freq(p, "marker")
#'
#' @export
calc_geno_freq <-
    function(probs, by=c("individual", "marker"), omit_x=TRUE)
{
    by <- match.arg(by)

    is_x_chr <- attr(probs, "is_x_chr")
    if(any(is_x_chr) && !all(is_x_chr) && omit_x) {
        probs <- probs[,!is_x_chr]
        is_x_chr <- is_x_chr[!is_x_chr]
        omit_x <- FALSE
    }

    if(any(is_x_chr) && any(!is_x_chr)) { # some autosome, some X chr
        return( list(A=calc_geno_freq(probs[,!is_x_chr], by, FALSE),
                     X=calc_geno_freq(probs[,is_x_chr], by, FALSE)) )
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

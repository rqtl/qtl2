#' Join genome scan results for different phenotypes.
#'
#' Join multiple \code{\link{scan1}} results for different phenotypes;
#' must have the same map.
#'
#' @param ... Genome scan objects as produced by \code{\link{scan1}}.
#' Must have the same map.
#'
#' @return A single genome scan object with the lod score columns
#' combined as multiple columns.
#'
#' @details If components \code{addcovar}, \code{Xcovar},
#' \code{intcovar}, \code{weights} do not match between objects, we
#' omit this information.
#'
#' If \code{hsq} present but has differing numbers of rows, we omit this information.
#'
#' If \code{snpinfo} present, they much match.
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' phe1 <- grav2$pheno[,1,drop=FALSE]
#' phe2 <- grav2$pheno[,2,drop=FALSE]
#'
#' out1 <- scan1(probs, phe1) # phenotype 1
#' out2 <- scan1(probs, phe2) # phenotype 2
#' out <- cbind(out1, out2)
#'
#' @export
cbind.scan1 <-
    function(...)
{
    args <- list(...)

    # to cbind: lod, coef, SE
    # to c(): sample_size
    # check if matches: map, snpinfo
    # drop if not matching: addcovar, Xcovar, intcovar, weights

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that maps match
    must_match <- c("map", "snpinfo", "is_x_chr")
    for(i in 2:length(args)) {
        for(obj in must_match) {
            if(!is_same(result[[obj]], args[[i]][[obj]]))
                stop("Input objects 1 and ", i, " have different ", obj)
        }
    }

    # cbind the lod (don't consider coef/SE for this, and require that lod are present)
    for(i in 2:length(args)) {
        if(!("lod" %in% names(args[[1]])) || !("lod" %in% names(args[[i]])))
            stop("Not all inputs have LOD score columns")
        if(!is_same(rownames(result$lod), rownames(args[[i]]$lod)))
            stop("Input objects 1 and ", i, " have different lod rownames.")
        result$lod <- cbind(result$lod, args[[i]]$lod)
    }

    # combine sample_size
    for(i in 2:length(args))
        result$sample_size <- c(result$sample_size, args[[i]]$sample_size)

    other_stuff <- c("addcovar", "Xcovar", "intcovar", "weights")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(result[[obj]], args[[i]][[obj]]))
                result[[obj]] <- NULL
        }
    }

    # hsq must have same rows or is omitted
    for(i in 2:length(args)) {
        if(!("hsq" %in% names(args[[1]])) && !("hsq" %in% names(args[[i]]))) next # not present
        if(!("hsq" %in% names(args[[1]])) || !("hsq" %in% names(args[[i]])))
            result$hsq <- NULL
        if(!is_same(rownames(result$hsq), rownames(args[[i]]$hsq)))
            result$hsq <- NULL
        result$hsq <- cbind(result$hsq, args[[i]]$hsq)
    }

    result
}

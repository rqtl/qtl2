#' Join genome scan results for different chromosomes.
#'
#' Join multiple \code{\link{scan1}} results for different
#' chromosomes; must have the same set of lod score column.
#'
#' @param ... Genome scan objects as produced by \code{\link{scan1}}.
#' Must have the same lod score columns.
#'
#' @return A single genome scan object with the results for
#' different sets of chromosomes combined.
#'
#' @details If components \code{addcovar}, \code{Xcovar},
#' \code{intcovar}, \code{weights}, \code{sample_size} do not match
#' between objects, we omit this information.
#'
#' If \code{hsq} present, we simply \code{rbind()} the contents.
#'
#' @examples
#' library(qtl2geno)
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
#' phe <- grav2$pheno[,1,drop=FALSE]
#'
#' out1 <- scan1(probs[,1], phe) # chr 1
#' out2 <- scan1(probs[,5], phe) # chr 5
#' out <- rbind(out1, out2)
#'
#' @export
rbind.scan1 <-
    function(...)
{
    args <- list(...)

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # to rbind: lod, coef, SE, hsq
    # to c(): map, snpinfo, is_x_chr
    # drop if not matching: addcovar, Xcovar, intcovar, weights, sample_size

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # cbind the lod, coef, SE
    main_stuff <- c("lod", "coef", "SE", "hsq")
    for(i in 2:length(args)) {
        for(obj in main_stuff) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not present in all inputs")
            if(!is_same(colnames(result[[obj]]), colnames(args[[i]][[obj]])))
                stop("Input objects 1 and ", i, " have different ", obj, " columns.")
            result[[obj]] <- rbind(result[[obj]], args[[i]][[obj]])
        }
    }

    # combine map, snpinfo, is_x_chr
    to_combine <- c("map", "snpinfo", "is_x_chr")
    for(i in 2:length(args)) {
        for(obj in to_combine) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not present in all inputs")
            result[[obj]] <- c(result[[obj]], args[[i]][[obj]])
        }
    }

    other_stuff <- c("addcovar", "Xcovar", "intcovar", "weights", "sample_size")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(result[[obj]], args[[i]][[obj]]))
                result[[obj]] <- NULL
        }
    }

    result
}

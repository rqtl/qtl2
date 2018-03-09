#' Join genome scan results for different chromosomes.
#'
#' Join multiple [scan1()] results for different
#' chromosomes; must have the same set of lod score column.
#'
#' @md
#'
#' @param ... Genome scan objects as produced by [scan1()].
#' Must have the same lod score columns.
#'
#' @return A single genome scan object with the results for
#' different sets of chromosomes combined.
#'
#' @details If components `addcovar`, `Xcovar`,
#' `intcovar`, `weights`, `sample_size` do not match
#' between objects, we omit this information.
#'
#' If `hsq` present, we simply `rbind()` the contents.
#'
#' @seealso [cbind.scan1()], [scan1()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
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
    if(length(args) == 1) return(result)


    # to rbind: main data, attributes SE, hsq
    # drop if not matching: sample_size

    # grab attributes
    args_attr <- lapply(args, attributes)

    result <- unclass(args[[1]])

    # cbind the main data
    for(i in 2:length(args)) {
        if(!is_same(colnames(result), colnames(args[[i]])))
           stop("Input objects 1 and ", i, " have different columns.")
           result <- rbind(result, unclass(args[[i]]))
    }

    # rbind attributes
    to_rbind <- c("hsq", "SE")
    for(i in 2:length(args)) {
        for(obj in to_rbind) {
            if(is.null(args_attr[[1]][[obj]]) && is.null(args_attr[[i]][[obj]])) next # not present
            if(is.null(args_attr[[1]][[obj]]) || is.null(args_attr[[i]][[obj]]))
                stop(obj, " not present in all inputs")
            args_attr[[1]][[obj]] <- rbind(args_attr[[1]][[obj]], args_attr[[i]][[obj]])
        }
    }

    # drop if not matching
    other_stuff <- c("sample_size")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(args_attr[[1]][[obj]], args_attr[[i]][[obj]]))
                args_attr[[1]][[obj]] <- NULL
        }
    }

    # copy attributes
    for(obj in c("sample_size", "hsq", "SE"))
        attr(result, obj) <- args_attr[[1]][[obj]]
    class(result) <- class(args[[1]])

    result
}

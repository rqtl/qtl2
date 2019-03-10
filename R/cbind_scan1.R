#' Join genome scan results for different phenotypes.
#'
#' Join multiple [scan1()] results for different phenotypes;
#' must have the same map.
#'
#' @param ... Genome scan objects as produced by [scan1()].
#' Must have the same map.
#'
#' @return A single genome scan object with the lod score columns
#' combined as multiple columns.
#'
#' @details If components `addcovar()`, `Xcovar`,
#' `intcovar`, `weights` do not match between objects, we
#' omit this information.
#'
#' If `hsq` present but has differing numbers of rows, we omit this information.
#'
#' @seealso [rbind.scan1()], [scan1()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
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
    if(length(args) == 1) return(result)

    # to cbind: main data, attributes SE, hsq
    # to c(): attribute sample_size

    # grab attributes
    args_attr <- lapply(args, attributes)

    result <- unclass(args[[1]])

    # cbind the data
    for(i in 2:length(args)) {
        if(!is_same(rownames(result), rownames(args[[i]])))
            stop("Input objects 1 and ", i, " have different rownames.")
        result <- cbind(result, unclass(args[[i]]))
    }

    # combine sample_size
    for(obj in c("sample_size")) {
        for(i in 2:length(args)) {
            drop_arg <- FALSE
            if(is.null(args_attr[[1]][[obj]]) && is.null(args_attr[[i]][[obj]])) next
            if(is.null(args_attr[[1]][[obj]]) || is.null(args_attr[[i]][[obj]])) {
                drop_arg <- TRUE
            } else {
                args_attr[[1]][[obj]] <- c(args_attr[[1]][[obj]], args_attr[[i]][[obj]])
            }
        }
        if(drop_arg) {
            warning('attribute "', obj, '" not present in all objects')
            args_attr[[1]][[obj]] <- NULL
        }

    }

    # attributes to combine by cbind()
    drop_arg <- c(hsq=FALSE, SE=FALSE)
    diff_rn_arg <- c(hsq=FALSE, SE=FALSE)
    for(obj in c("hsq", "SE")) {
        drop_arg <- FALSE
        drop_rn_arg <- FALSE
        for(i in 2:length(args)) {
            if(is.null(args_attr[[1]][[obj]]) && is.null(args_attr[[i]][[obj]])) next
            if(is.null(args_attr[[1]][[obj]]) || is.null(args_attr[[i]][[obj]])) {
                drop_arg <- TRUE
            } else if(!is_same(rownames(args_attr[[1]][[obj]]), rownames(args_attr[[i]][[obj]]))) {
                diff_rn_arg <- TRUE
            } else {
                args_attr[[1]][[obj]] <- cbind(args_attr[[1]][[obj]], args_attr[[i]][[obj]])
            }
        }
        if(drop_arg) {
            warning('attribute "', obj, '" not present in all objects')
            args_attr[[1]][[obj]] <- NULL
        } else if(drop_rn_arg) {
            warning('attribute "', obj, '" does not have the same row names in all objects')
            args_attr[[1]][[obj]] <- NULL
        }
    }

    for(obj in c("sample_size", "hsq", "SE"))
        attr(result, obj) <- args_attr[[1]][[obj]]
    class(result) <- class(args[[1]])

    result
}

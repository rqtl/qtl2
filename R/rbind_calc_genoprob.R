#' Join genotype probabilities for different individuals
#'
#' Join multiple genotype probability objects, as produced by
#' \code{\link{calc_genoprob}} for different individuals.
#'
#' @param ... Genotype probability objects as produced by
#' \code{\link{calc_genoprob}}. Must have the same set of markers and
#' genotypes.
#'
#' @return A single genotype probability object.
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probsA <- calc_genoprob(grav2[1:5,], step=1, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[6:12,], step=1, error_prob=0.002)
#' probs <- rbind(probsA, probsB)
#'
#' @export
rbind.calc_genoprob <-
    function(...)
{
    args <- list(...)

    # to rbind: probs, sex, cross_info
    # to pass through (must match): map, grid, crosstype, is_x_chr, alleles, alleleprobs, step, off_end, stepwidth

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    nested_stuff <- c("map", "grid")
    other_stuff <- c("crosstype", "is_x_chr", "alleles", "alleleprobs",
                     "step", "off_end", "stepwidth")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(args[[1]][[obj]], args[[i]][[obj]]))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
        for(obj in nested_stuff) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!(obj %in% names(args[[1]])) || !(obj %in% names(args[[i]])))
                stop(obj, " not prsent in all inputs")
            if(!is_same(names(args[[1]][[obj]]), names(args[[i]][[obj]])))
                stop("Input objects 1 and ", i, " differ in their ", obj)
            for(chr in seq(along=args[[1]][[obj]])) {
                if(!is_same(args[[1]][[obj]][[chr]], args[[i]][[obj]][[chr]]))
                    stop("Input objects 1 and ", i, " differ in their ", obj,
                         " on chromosome ", chr)
            }
        }
    }

    # create space for result
    for(obj in c("probs", "draws")) {
        if(obj %in% names(args[[1]])) {
            nind <- vapply(args, function(a) nrow(a[[obj]][[1]]), 1)
            totind <- sum(nind)
            index <- split(1:totind, rep(seq(along=nind), nind))

            result[[obj]] <- vector("list", length(args[[1]][[obj]]))
            names(result[[obj]]) <- names(args[[1]][[obj]])
            for(chr in names(result[[obj]])) {
                dimn <- dimnames(args[[1]][[obj]][[chr]])
                dimv <- dim(args[[1]][[obj]][[chr]])
                result[[obj]][[chr]] <- array(dim=c(totind, dimv[2], dimv[3]))
                dimnames(result[[obj]][[chr]]) <- list(paste(1:totind), dimn[[2]], dimn[[3]])
            }
        }
    }


    # paste stuff together
    nested_stuff <- c("probs", "draws")
    other_stuff <- c("sex", "cross_info")
    for(i in 1:length(args)) {
        for(obj in c("probs", "draws")) {
            if(!(obj %in% names(args[[1]])) && !(obj %in% names(args[[i]]))) next # not present
            if(!is_same(names(args[[1]][[obj]]), names(args[[i]][[obj]])))
                stop("Input objects 1 and ", i, " differ in the their ", obj)
            for(chr in names(args[[1]][[obj]])) {
                dimn1 <- dimnames(args[[1]][[obj]][[chr]])
                dimni <- dimnames(args[[i]][[obj]][[chr]])
                if(!is_same(dimn1[-1], dimni[-1]))
                    stop("Input objects 1 and ", i, " differ in the their ", obj,
                         " on chromosome ", chr)

                result[[obj]][[chr]][index[[i]],,] <- args[[i]][[obj]][[chr]]
                rownames(result[[obj]][[chr]])[index[[i]]] <- rownames(args[[i]][[obj]][[chr]])
            }
        }
    }

    for(i in 2:length(args)) {
        if(!("sex" %in% names(result)) && !("sex" %in% names(args[[i]]))) next
        if(!("sex" %in% names(result) && "sex" %in% names(args[[i]])))
            stop("sex present in only some input objects")

        if(!("cross_info" %in% names(result)) && !("cross_info" %in% names(args[[i]]))) next
        if(!("cross_info" %in% names(result) && "cross_info" %in% names(args[[i]])))
            stop("cross_info present in only some input objects")
        if(!is_same(ncol(result$cross_info), ncol(args[[i]]$cross_info)))
            stop("input objects have varying numbers of cross_info columns")
        result$cross_info <- rbind(result$cross_info, args[[i]]$cross_info)
    }

    result
}

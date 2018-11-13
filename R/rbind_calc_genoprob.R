#' Join genotype probabilities for different individuals
#'
#' Join multiple genotype probability objects, as produced by
#' [calc_genoprob()], for the same set of markers and genotypes but for
#' different individuals.
#'
#' @md
#'
#' @param ... Genotype probability objects as produced by
#' [calc_genoprob()]. Must have the same set of markers and
#' genotypes.
#'
#' @return A single genotype probability object.
#'
#' @seealso [cbind.calc_genoprob()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probsA <- calc_genoprob(grav2[1:5,], map, error_prob=0.002)
#' probsB <- calc_genoprob(grav2[6:12,], map, error_prob=0.002)
#' probs <- rbind(probsA, probsB)
#'
#' @export
rbind.calc_genoprob <-
    function(...)
{
    args <- list(...)

    # to rbind: the data
    # to pass through (must match): crosstype, is_x_chr, alleles, alleleprobs

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    other_stuff <- c("crosstype", "is_x_chr", "alleles", "alleleprobs")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
    }

    # create space for result
    nind <- vapply(args, function(a) nrow(a[[1]]), 1)
    totind <- sum(nind)
    index <- split(1:totind, rep(seq(along=nind), nind))

    result <- vector("list", length(args[[1]]))
    names(result) <- names(args[[1]])
    for(chr in names(result)) {
        dimn <- dimnames(args[[1]][[chr]])
        dimv <- dim(args[[1]][[chr]])
        result[[chr]] <- array(dim=c(totind, dimv[2], dimv[3]))
        dimnames(result[[chr]]) <- list(paste(1:totind), dimn[[2]], dimn[[3]])
    }

    # paste stuff together
    for(i in 1:length(args)) {
        if(!is_same(names(args[[1]]), names(args[[i]])))
            stop("Input objects 1 and ", i, " have different chromosome names")
        for(chr in names(args[[1]])) {
            dimn1 <- dimnames(args[[1]][[chr]])
            dimni <- dimnames(args[[i]][[chr]])
            if(!is_same(dimn1[-1], dimni[-1]))
                stop("Input objects 1 and ", i, " differ in shape on chromosome ", chr)

            result[[chr]][index[[i]],,] <- args[[i]][[chr]]
            rownames(result[[chr]])[index[[i]]] <- rownames(args[[i]][[chr]])
        }
    }

    # paste in the attributes
    for(obj in other_stuff)
        attr(result, obj) <- attr(args[[1]], obj)
    class(result) <- class(args[[1]])

    result
}

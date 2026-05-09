#' Join phased geno results for different individuals
#'
#' Join phased genotype objects, as produced by [guess_phase()],
#' for the same set of markers but for different individuals.
#'
#' @param ... Imputed genotype objects as produced by
#' [guess_phase()]. Must have the same set of markers.
#'
#' @return An object of class `"phasedgeno"`, like the input; see [guess_phase()].
#'
#' @seealso [cbind.phasedgeno()], [subset.phasedgeno()], [guess_phase()]
#'
#' @examples
#' \dontrun{
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' m <- maxmarg(pr)
#' phA <- guess_phase(DOex[1:10,], m[1:10,])
#' phB <- guess_phase(DOex[11:20,], m[11:20,])
#' ph <- rbind(phA, phB)
#' }
#'
#' @export
rbind.phasedgeno <-
    function(...)
{
    args <- list(...)

    # to rbind: the data itself
    # to pass through (must match): crosstype, is_x_chr, alleles

    result <- args[[1]]
    if(length(args) == 1) return(result)

    # check that things match
    other_stuff <- c("crosstype", "is_x_chr", "alleles")
    for(i in 2:length(args)) {
        for(obj in other_stuff) {
            if(!is_same(attr(args[[1]], obj), attr(args[[i]], obj)))
                stop("Input objects 1 and ", i, " differ in their ", obj)
        }
    }

    for(i in 2:length(args)) {
        if(length(result) != length(args[[i]]) ||
           !all(names(result) == names(args[[i]])))
            stop("Input arguments have different chromosomes")

        for(s in seq(along=result)) {
            if(!is_same(ncol(result[[s]]), ncol(args[[i]][[s]])))
                stop("input objects have varying numbers of geno columns on chr ", names(result)[s])
            if(!is_same(dim(result[[s]])[3], dim(args[[i]][[s]])[3]))
                stop("input objects have varying numbers of haplotype columns on chr ", names(result)[s])
            temp <- array(dim=c(nrow(result[[s]])+nrow(args[[i]][[s]]),
                                dim(result[[s]])[2:3]))
            dimnames(temp) <- list(c(rownames(result[[s]]), rownames(args[[i]][[s]])),
                                   colnames(result[[s]]), dimnames(result[[s]])[[3]])
            for(j in 1:dim(result[[s]])[3]) temp[,,j] <- rbind(result[[s]][,,j], args[[i]][[s]][,,j])
            result[[s]] <- temp
        }
    }

    result
}

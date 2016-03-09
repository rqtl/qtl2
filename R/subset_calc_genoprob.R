# subset calc_genoprob objects

#' Subsetting genotype probabilities
#'
#' Pull out a specified set of individuals and/or chromosomes from
#' the results of \code{\link{calc_genoprob}}.
#'
#' @param x Genotype probabilities as output from \code{\link{calc_genoprob}}.
#' @param ind A vector of individuals: numeric indices, logical
#' values, or character string IDs
#' @param chr A vector of chromosomes: logical values, or character
#' string IDs. Numbers are interpreted as character string IDs.
#' @param ... Ignored.
#'
#' @return The input genotype probabilities, with the selected
#' individuals and/or chromsomes.
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' \dontshow{grav2 <- grav2[1:8,c(1,2)]}
#' pr <- calc_genoprob(grav2)
#' # keep just individuals 1:5, chromosome 2
#' prsub <- pr[1:5,2]
#' # keep just chromosome 2
#' prsub2 <- pr[,2]
subset.calc_genoprob <-
    function(x, ind, chr, ...)
{
    if((missing(ind) || is.null(ind)) &&
       (missing(chr) || is.null(chr)))
        stop("You must specify either ind or chr.")

    chrID <- names(x$map)
    n_chr <- length(chrID)
    if(!missing(chr) && !is.null(chr)) {
        if(is.logical(chr)) {
            if(length(chr) != n_chr) {
                stop("length(chr) [", length(chr), "] != no. chr in x [",
                     n_chr, "]")
            }
            chr <- chrID[chr] # convert to character strings
        }
        else {
            chrindex <- match(chr, chrID)
            if(any(is.na(chrindex))) {
                stop("Some chr not found: ",
                     paste(chr[is.na(chrindex)], collapse=", "))
            }
            chr <- chrID[chrindex] # character strings
        }
        if(length(chr) == 0)
            stop("Must retain at least one chromosome.")

        to_sub <- c("probs", "draws", "map", "is_x_chr", "grid", "snpinfo") # draws is here, to also deal with sim_geno objects
        for(a in to_sub) {
            if(a %in% names(x)) {
                x[[a]] <- x[[a]][chr]
            }
        }
    }

    if("probs" %in% names(x))
        indID <- rownames(x$probs[[1]])
    else if("draws" %in% names(x))
        indID <- rownames(x$draws[[1]])
    else stop("Neither probs no draws found.")

    n_ind <- length(indID)
    if(!missing(ind) && !is.null(ind)) {
        if(is.logical(ind)) {
            if(length(ind) != n_ind) {
                stop("length(ind) [", length(ind), "] != no. ind in x [",
                     n_ind, "]")
            }
            ind <- indID[ind] # convert to character strings
        }
        else if(is.numeric(ind)) {
            if(any(ind < 1 || ind > n_ind))
                stop("Numeric ind out of allowed range [1 - ", n_ind, "]")
            ind <- indID[ind] # convert to character strings
        }
        else {
            indindex <- match(ind, indID)
            if(any(is.na(indindex))) {
                stop("Some individuals not found: ",
                     paste(ind[is.na(indindex)], collapse=", "))
            }
            ind <- indID[indindex]
        }
        if(length(ind) == 0)
            stop("Must retain at least one individual.")

        for(a in c("probs", "draws")) {
            if(a %in% names(x)) {
                for(i in names(x$map)) # loop over chromosomes
                    x[[a]][[i]] <- x[[a]][[i]][ind,,,drop=FALSE]
            }
        }
        if(!is.null(x$cross_info))
            x$cross_info <- x$cross_info[ind, , drop=FALSE]
    }

    x
}

#' @export
#' @rdname subset.calc_genoprob
`[.calc_genoprob` <-
    function(x, ind, chr)
    subset(x, ind, chr)

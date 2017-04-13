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
    function(x, ind=NULL, chr=NULL, ...)
{
    if(is.null(ind) && is.null(chr))
        stop("You must specify either ind or chr.")

    if(!is.null(chr)) {
        chr <- subset_chr(chr, names(x))

        if(length(chr) == 0)
            stop("Must retain at least one chromosome.")

        cl <- class(x)
        class(x) <- "list"

        attr_to_sub <- c("is_x_chr")
        attr_to_keep <- c("crosstype", "alleles", "alleleprobs")
        x_attr <- attributes(x)
        x_attrnam <- names(x_attr)
        x <- x[chr]
        for(a in attr_to_sub) {
            if(a %in% x_attrnam) {
                attr(x, a) <- x_attr[[a]][chr]
            }
        }
        for(a in attr_to_keep) {
            if(a %in% attr_to_keep)
                attr(x, a) <- x_attr[[a]]
        }
        class(x) <- cl
    }

    if(!is.null(ind)) {
        if("calc_genoprob" %in% class(x))
            all_ind <- dimnames(x)[[1]]
        else all_ind <- rownames(x[[1]])

        ind <- subset_ind(ind, all_ind)

        if(length(ind) == 0)
            stop("Must retain at least one individual.")

        cl <- class(x)
        class(x) <- "list"

        for(i in names(x)) # loop over chromosomes
            x[[i]] <- x[[i]][ind,,,drop=FALSE]

        class(x) <- cl
    }

    x
}

#' @export
#' @rdname subset.calc_genoprob
`[.calc_genoprob` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)

# subset calc_genoprob objects
# ...need this in case qtl2geno isn't available [eye roll]

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

`[.calc_genoprob` <-
    function(x, ind=NULL, chr=NULL)
    subset(x, ind, chr)

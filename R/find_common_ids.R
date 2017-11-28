# match ids
#
# takes two vectors of character strings
# returns a list with
#   first,second = indexes in both vectors, to make them aligned
#   firstonly = indexes in first vector that are not in the second
#   secondonly = indexes in second vector that are not in the first
match_ids <-
    function(x,y)
{
    mxy <- match(x, y)
    myx <- match(y, x)

    result <- list(first=which(!is.na(mxy)),
                   second=mxy[!is.na(mxy)],
                   firstonly=which(is.na(mxy)),
                   secondonly=which(is.na(myx)))

    if(length(result$first) > 0)
        names(result$first) <- x[result$first]
    if(length(result$firstonly) > 0)
        names(result$firstonly) <- x[result$firstonly]
    if(length(result$second) > 0)
        names(result$second) <- y[result$second]
    if(length(result$secondonly) > 0)
        names(result$secondonly) <- y[result$secondonly]

    result
}

# find common ids across a set of vectors of character string
# return value: ids that are present in *all* of the inputs
find_common_ids <-
    function(...)
{
    inp <- list(...)
    n_inp <- length(inp)

    if(n_inp==0)
        stop("No input")

    if(n_inp==1) {
        if(length(inp[[1]])==0) return(character(0))
        return(inp[[1]])
    }

    reduced <- inp[[1]]
    for(i in 2:n_inp) {
        m <- match_ids(reduced, inp[[i]])
        reduced <- reduced[m$first]
    }

    if(length(reduced)==0) return(character(0))

    reduced
}

# align genotypes, is_female and cross_info
# (because we'll do this alot)
align_geno_sex_cross <-
    function(geno, is_female, cross_info, hold_off_warning=FALSE)
{
    if(is.list(geno)) {
        give_warning <- FALSE
        dropped_ind <- NULL
        for(i in seq(along=geno)) {
            orig_nind <- nrow(geno[[i]])
            orig_ind <- rownames(geno[[i]])
            tmp <- align_geno_sex_cross(geno[[i]],
                                        is_female,
                                        cross_info,
                                        hold_off_warning=TRUE)
            if(orig_nind > nrow(tmp$geno)) {
                new_ind <- rownames(tmp$geno)
                dropped_ind <- unique(c(dropped_ind,
                                        orig_ind[!(orig_ind %in% new_ind)]))
                give_warning <- TRUE
            }
            is_female <- tmp$is_female
            cross_info <- tmp$cross_info
            geno[[i]] <- tmp$geno
        }
        if(give_warning)
            warning("Omitting genotypes for ", length(dropped_ind),
                    ifelse(length(dropped_ind)==1, " individual", " individuals"),
                    " with no sex/cross info.")

        return(list(geno=geno, is_female=is_female, cross_info=cross_info))
    }

    ind <- rownames(geno)
    is_female <- handle_null_isfemale(is_female, ind)
    cross_info <- handle_null_crossinfo(cross_info, ind)

    keep <- find_common_ids(rownames(geno), names(is_female), rownames(cross_info))

    if(length(keep) < length(ind) && !hold_off_warning)
        warning("Omitting genotypes for ", length(ind) - length(keep),
                ifelse(length(ind)-length(keep)==1, " individual", " individuals"),
                " with no sex/cross info.")

    list(geno=geno[keep,,drop=FALSE],
         is_female=is_female[keep],
         cross_info=cross_info[keep,,drop=FALSE])
}

# subset kinship matrix by individuals or chromosomes
subset_kinship <-
    function(kinship, ind=NULL, chr=NULL)
{
    if(is.null(kinship)) return(NULL)

    # already decomposed
    if(is_kinship_decomposed(kinship)) {
        # maybe it's not really needing to be subset by individual
        if(is_kinship_list(kinship)) {
            k_ind <- rownames(kinship$vectors[[1]])
        } else {
            k_ind <- rownames(kinship$vectors)
        }
        new_ind <- subset_ind(ind, k_ind)
        if(length(new_ind) != length(k_ind) ||
           any(new_ind != k_ind)) {
            stop("Can't subset decomposed kinship matrices by individual")
        }
        # otherwise, not really being subset so can ignore ind argument

        if(!is_kinship_list(kinship)) { # can ignore chr argument
            return(kinship)
        }

        if(!is.null(chr)) {
            chr <- subset_chr(chr, names(kinship))
            kinship <- kinship[chr]
        }
        if(length(kinship)==1) return(kinship[[1]])
        return(kinship)
    }

    # not yet decomposed
    if(!is.null(chr) && is.list(kinship)) {
        chr <- subset_chr(chr, names(kinship))
        kinship <- kinship[chr]

        if(is.list(kinship) && length(kinship) == 1)
            kinship <- kinship[[1]]
    }

    if(!is.null(ind)) {
        if(is.list(kinship)) {
            ind <- subset_ind(ind, rownames(kinship[[1]]))
            kinship <- lapply(kinship, function(a) a[ind, ind])
        }
        else {
            ind <- subset_ind(ind, rownames(kinship))
            kinship <- kinship[ind,ind]
        }
    }

    kinship
}

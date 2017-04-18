# subset kinship matrix by individuals or chromosomes
subset_kinship <-
    function(kinship, ind=NULL, chr=NULL)
{
    if(is.null(kinship)) return(NULL)

    # already decomposed
    decomped <- attr(kinship, "eigen_decomp")
    if(!is.null(decomped) && decomped) {
        if(!is.null(ind))
            stop("Can't subset decomposed kinship matrices by individual")

        if(is.list(kinship) && length(kinship)==2 &&
           all(names(kinship) == c("values", "vectors")))
            return(kinship)


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

# function for subsetting by individual
#
# ind = logical, numeric, or character strings
# all_ind = vector of character strings with all individual names
#
# result = vector of character strings (the selected ones)
subset_ind <-
    function(ind, all_ind, allow_logical=TRUE)
{
    if(is.logical(ind)) {
        if(!allow_logical) stop("ind can't be logical if different individuals in geno, pheno, covar")

        if(length(ind) != length(all_ind))
            stop("ind is logical but length [", length(ind), "] != no. ind in cross [",
                 length(all_ind), "]")
        ind <- all_ind[ind]
    }
    if(is.numeric(ind)) { # treat as numeric indexes
        if(!allow_logical) stop("ind can't be numeric if different individuals in geno, pheno, covar")
        if(any(ind < 0)) { # deal with negatives
            if(!all(ind < 0)) stop("Can't mix negative and positive ind subscripts")
                ind <- (seq_along(all_ind))[ind]
        }

        if(any(ind < 1 | ind > length(all_ind))) {
            ind <- ind[ind >= 1 & ind <= length(all_ind)]
            if(length(ind)==0) stop("All ind out of range")
            warning("some ind out of range [1,", length(all_ind), "]")
        }

        ind <- all_ind[ind]
    }

    # now treat ind as character strings
    ind <- as.character(ind)

    # look for negatives; turn to positives
    if(!any(grepl("^\\-", all_ind)) && any(grepl("^\\-", ind))) { # if "-" used in actual IDs, don't allow negative subscripts
        if(!all(grepl("^\\-", ind)))
            stop("Can't mix negative and positive ind subscripts")
        ind <- sub("^\\-", "", ind)
        if(!all(ind %in% all_ind)) {
            if(!any(ind %in% all_ind)) stop("None of the individuals in cross")
            warning("Some individuals not in cross: ", paste(ind[!(ind %in% all_ind)], collapse=", "))
        }
        ind <- all_ind[!(all_ind %in% ind)]
    }

    if(!all(ind %in% all_ind)) {
        if(!any(ind %in% all_ind)) stop("None of the individuals in cross")
        warning("Some ind not in cross: ", paste(ind[!(ind %in% all_ind)], collapse=", "))
        ind <- ind[ind %in% ind]
    }

    ind
}

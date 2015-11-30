# functions to deal with covariates

# force interactive covariates into the additive covariate matrix
#
# addcovar and intcovar are two matrices
# intcovar columns should all be within the addcovar columns
force_intcovar <-
    function(addcovar=NULL, intcovar=NULL)
{
    if(is.null(intcovar)) # no intcovar, so return addcovar w/o change
        return(addcovar)

    if(is.null(addcovar)) # no addcovar, so just return intcovar
        return(intcovar)

    if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
    if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)

    stopifnot(nrow(addcovar) == nrow(intcovar))
    n.addcovar <- ncol(addcovar)
    n.intcovar <- ncol(intcovar)

    # look for matching columns
    full <- cbind(addcovar, intcovar)
    has_match <- find_matching_cols(full)
    if(any(has_match > 0))
        full <- full[, has_match<0, drop=FALSE]

    full
}

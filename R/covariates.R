# functions to deal with covariates

# force interactive covariates into the additive covariate matrix
#
# addcovar and intcovar are two matrices
# intcovar columns should all be within the addcovar columns
# tol is tolerance for determining matching columns
force_intcovar <-
    function(addcovar=NULL, intcovar=NULL, tol=1e-12)
{
    if(is.null(intcovar)) # no intcovar, so return addcovar w/o change
        return(addcovar)

    if(is.null(addcovar)) # no addcovar, so just return intcovar
        return(intcovar)

    if(!is.matrix(addcovar)) addcovar <- as.matrix(addcovar)
    if(!is.matrix(intcovar)) intcovar <- as.matrix(intcovar)

    # IDs in both; omitting any individuals with missing values
    ids <- get_common_ids(addcovar[complete.cases(addcovar),,drop=FALSE],
                          intcovar[complete.cases(intcovar),,drop=FALSE])

    # look for matching columns, having reduced to common individuals
    full <- cbind(addcovar[ids,,drop=FALSE], intcovar[ids,,drop=FALSE])
    has_match <- find_matching_cols(full, tol)
    if(any(has_match > 0))
        full <- full[, has_match<0, drop=FALSE]

    if(ncol(full)==0) return(NULL)

    full
}

# drop linearly dependent columns
# if intercept=TRUE, add intercept before checking and then remove afterwards
drop_depcols <-
    function(covar=NULL, add_intercept=FALSE, tol=1e-12)
{
    if(is.null(covar)) return(covar)

    if(!is.matrix(covar)) covar <- as.matrix(covar)
    if(add_intercept) covar <- cbind(rep(1, nrow(covar)), covar)

    if(ncol(covar) <= 1) return(covar)

    # deal with NAs by omitting those rows before
    depcols <- sort(find_lin_indep_cols(covar[complete.cases(covar),,drop=FALSE], tol))

    if(add_intercept) {
        # FIX_ME: assuming here that intercept (first column) will always be included
        depcols <- depcols[-1]
    }
    if(length(depcols)==0) return(NULL)

    covar[, depcols, drop=FALSE]
}

# drop columns from X covariates that are already in addcovar
drop_xcovar <-
    function(covar=NULL, Xcovar=NULL, tol=1e-12)
{
    if(is.null(Xcovar) || is.null(covar)) return(Xcovar)

    if(!is.matrix(covar)) covar <- as.matrix(covar)
    if(!is.matrix(Xcovar)) Xcovar <- as.matrix(Xcovar)

    # IDs in both; omitting any individuals with missing values
    ids <- get_common_ids(covar[complete.cases(covar),,drop=FALSE],
                          Xcovar[complete.cases(Xcovar),,drop=FALSE])

    # look for matching columns, having reduced to common individuals
    matches <- find_matching_cols(cbind(covar[ids,], Xcovar[ids,]), tol)[-(1:ncol(covar))]

    if(all(matches > 0)) return(NULL)

    # drop the columns with matches
    Xcovar[,matches<0,drop=FALSE]
}

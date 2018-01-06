# functions to deal with covariates

# force interactive covariates into the additive covariate matrix
#
# addcovar and intcovar are two matrices
# intcovar columns should all be within the addcovar columns
# tol is tolerance for determining matching columns
#' @importFrom stats complete.cases
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
#' @importFrom stats complete.cases
drop_depcols <-
    function(covar=NULL, add_intercept=FALSE, tol=1e-12)
{
    if(is.null(covar)) return(covar)

    if(!is.matrix(covar)) covar <- as.matrix(covar)
    if(add_intercept) covar <- cbind(rep(1, nrow(covar)), covar)

    if(ncol(covar) <= 1) return(covar)

    X <- covar[complete.cases(covar),,drop=FALSE]

    # deal with NAs by omitting those rows before
    indep_cols <- sort(find_lin_indep_cols(X, tol))

    if(add_intercept) {

        target_ncol <- length(indep_cols)

        while(!(1 %in% indep_cols)) {
            # don't want to omit the first column (the intercept)
            # need to work harder...
            #  - drop one column at a time other the intercept
            #  - when you find a column that doesn't reduce the target number of columns, drop it
            #  - check again if the intercept is being retained; if not, repeat

            for(i in 2:ncol(X)) {
                indep_cols <- find_lin_indep_cols(X[,-i,drop=FALSE], tol)
                if(length(indep_cols) == target_ncol) {
                    X <- X[,-i,drop=FALSE]
                    break
                }
            }

            indep_cols <- sort(find_lin_indep_cols(X, tol))
        }

        # now drop the intercept
        indep_cols <- indep_cols[-1]
    }
    if(length(indep_cols)==0) return(NULL)

    covar[, indep_cols, drop=FALSE]
}

# drop columns from X covariates that are already in addcovar
#' @importFrom stats complete.cases
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
    matches <- find_matching_cols(cbind(covar[ids,], Xcovar[ids,]), tol)[-seq_len(ncol(covar))]

    if(all(matches > 0)) return(NULL)

    # drop the columns with matches
    Xcovar[,matches<0,drop=FALSE]
}

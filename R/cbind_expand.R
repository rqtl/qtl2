#' Combine matrices by columns, expanding and aligning rows
#'
#' This is like [base::cbind()] but using row names to align the rows and expanding
#' with missing values if there are rows in some matrices but not others.
#'
#' @md
#'
#' @param ... A set of matrices or data frames
#'
#' @return The matrices combined by columns, using row names to align the rows, and expanding with missing values if there are rows in some matrices but not others.
#'
#' @export
#'
#' @examples
#' df1 <- data.frame(x=c(1,2,3,NA,4), y=c(5,8,9,10,11), row.names=c("A", "B", "C", "D", "E"))
#' df2 <- data.frame(w=c(7,8,0,9,10), z=c(6,NA,NA,9,10), row.names=c("A", "B", "F", "C", "D"))
#' cbind_expand(df1, df2)
cbind_expand <-
    function(...)
{
    input <- list(...)

    if(length(input)<2) return(input[[1]])

    # check IDs
    id <- lapply(input, rownames)
    if(any(vapply(id, is.null, TRUE)))
        stop("All input matrices must have rownames (containing individual IDs)")
    if(any(vapply(id, function(a) length(unique(a)) != length(a), TRUE)))
        stop("Some input matrices have duplicate rownames")
    uid <- unique(unlist(id))

    for(i in seq(along=input)) {
        not_in <- !(uid %in% id[[i]])
        if(!any(not_in)) next
        missing <- uid[not_in]

        # create new rows
        new_rows <- matrix(NA, ncol=ncol(input[[i]]), nrow=length(missing))
        dimnames(new_rows) <- list(missing, colnames(input[[i]]))

        input[[i]] <- rbind(input[[i]], new_rows)
    }

    for(i in seq(along=input)) {
        if(i==1)
            input[[1]] <- input[[1]][uid,,drop=FALSE]
        else
            input[[1]] <- cbind(input[[1]], input[[i]][uid,,drop=FALSE])
    }

    input[[1]]
}


# this is the version that was used in qtl2scan, with no attempt to align the rows
cbind_expand_noalign <-
    function(...)
{
    input <- list(...)

    if(length(input)<2) return(input[[1]])

    if(!all(vapply(input, is.matrix, TRUE)))
        stop("Not all inputs are matrices")

    nrow <- vapply(input, nrow, 1)
    max_nrow <- max(nrow)

    for(i in seq_along(input)) {
        rownames(input[[i]]) <- NULL # strip off the rownames
        if(nrow[i] < max_nrow) { # need to pad
            new_rows <- matrix(NA, ncol=ncol(input[[i]]), nrow=max_nrow - nrow[i])
            colnames(new_rows) <- colnames(input[[i]])
            input[[i]] <- rbind(input[[i]], new_rows)
        }
    }

    do.call("cbind", input)
}

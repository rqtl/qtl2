#' Align objects by their row names
#'
#' Align a set of objects by their row names, omitting rows that are
#' not present in all objects.
#'
#' @param ... A set of objects: vectors, lists, matrices, data frames, and/or arrays
#' @param return_index If true, just return a vector with the
#' individuals that are in common across all objects.
#'
#' @return If \code{return_index=FALSE} (the default), a list
#' containing the objects reduced to the individuals in common across
#' all of them. If \code{return_index=TRUE}, just return a vector of
#' character strings for the individuals that are in common.
#'
#' @export
align_byrow <-
    function(..., return_index=FALSE)
{
    args <- list(...)
    if(length(args)==0) {
        if(return_index) return(character(0))
        else return(list())
    }

    # find the IDs in common across all
    id <- NULL
    for(i in seq(along=args)) {
        if(is.matrix(args[[i]]) || is.data.frame(args[[i]]) || is.array(args[[i]])) {
            if(length(dim(args[[i]])) > 3)
                stop("Can't handle arrays with >3 dimensions")
            these <- rownames(args[[i]])
        }
        else if(is.vector(args[[i]])) {
            these <- names(args[[i]])
        }
        else {
            stop("Not sure what to do with object of class ", class(args[[i]]))
        }

        if(length(unique(these)) != length(these))
            stop("Duplicate names in argument ", i)

        if(is.null(id)) id <- these
        else id <- id[id %in% these]
    }

    if(return_index) return(id)

    # subset the objects, preserving shape
    for(i in seq(along=args)) {
        if(is.matrix(args[[i]]) || is.data.frame(args[[i]]) || is.array(args[[i]])) {
            d <- length(dim(args[[i]]))
            if(d==2)
                args[[i]] <- args[[i]][id,,drop=FALSE]
            else if(d==3)
                args[[i]] <- args[[i]][id,,,drop=FALSE]
            else
                stop("Can't handle arrays with ", d, " dimensions")
        }
        else if(is.vector(args[[i]])) {
            args[[i]] <- args[[i]][id]
        }
        else {
            stop("Not sure what to do with object of class ", class(args[[i]]))
        }
    }

    args
}

#' Get common set of IDs from objects
#'
#' For a set objects with IDs as row names (or, for a vector, just
#' names), find the IDs that are present in all of the objects.
#'
#' @param ... A set of objects: vectors, lists, matrices, data frames,
#' and/or arrays. If one is a character vector with no names
#' attribute, it's taken to be a set of IDs, itself.
#'
#' @return A vector of character strings for the individuals that are
#' in common.
#'
#' @export
get_common_ids <-
    function(...)
{
    args <- list(...)
    if(length(args)==0) {
        return(character(0))
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
            if(is.character(args[[i]]) && is.null(names(args[[i]])))
                these <- args[[i]]
            else
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

    id
}

#' Get common set of IDs from objects
#'
#' For a set objects with IDs as row names (or, for a vector, just
#' names), find the IDs that are present in all of the objects.
#'
#' @param ... A set of objects: vectors, lists, matrices, data frames,
#' and/or arrays. If one is a character vector with no names
#' attribute, it's taken to be a set of IDs, itself.
#' @param complete.cases If TRUE, look at matrices and non-character
#' vectors and keep only individuals with no missing values.
#'
#' @return A vector of character strings for the individuals that are
#' in common.
#'
#' @details This is used (mostly internally) to align phenotypes,
#' genotype probabilities, and covariates in preparation for a genome
#' scan. The \code{complete.cases} argument is used to omit
#' individuals with any missing covariate values.
#'
#' @examples
#' x <- matrix(0, nrow=10, ncol=5); rownames(x) <- LETTERS[1:10]
#' y <- matrix(0, nrow=5, ncol=5);  rownames(y) <- LETTERS[(1:5)+7]
#' z <- LETTERS[5:15]
#' get_common_ids(x, y, z)
#'
#' x[8,1] <- NA
#' get_common_ids(x, y, z)
#' get_common_ids(x, y, z, complete.cases=TRUE)
#'
#' @export
get_common_ids <-
    function(..., complete.cases=FALSE)
{
    args <- list(...)
    if(length(args)==0) {
        return(character(0))
    }

    # find the IDs in common across all
    id <- NULL
    for(i in seq_along(args)) {
        if(is.null(args[[i]])) next

        if(is.matrix(args[[i]]) || is.data.frame(args[[i]]) || is.array(args[[i]])) {
            if(length(dim(args[[i]])) > 3)
                stop("Can't handle arrays with >3 dimensions")
            these <- rownames(args[[i]])
            if(complete.cases && (is.matrix(args[[i]]) || is.data.frame(args[[i]])))
                these <- these[complete.cases(args[[i]])]
        }
        else if(is.list(args[[i]]) && !is.null(rownames(args[[i]][[1]]))) {
            these <- rownames(args[[i]][[1]])
        }
        else if(is.vector(args[[i]])) {
            if(is.character(args[[i]]) && is.null(names(args[[i]])))
                these <- args[[i]]
            else {
                these <- names(args[[i]])
                if(complete.cases)
                    these <- these[!is.na(args[[i]])]
            }
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

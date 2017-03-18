#' Calculate chromosome lengths
#'
#' Calculate chromosome lengths for a map object
#'
#' @param map A list of vectors, each specifying locations of the markers.
#' @param collapse_to_AX If TRUE, collapse to the total lengths of the
#' autosomes and X chromosome.
#'
#' @return A vector of chromosome lengths. If
#' \code{collapse_to_AX=TRUE}, the result is a vector of length 2
#' (autosomal and X chromosome lengths).
#'
#' @details We take \code{diff(range(v))} for each vector, \code{v}.
#'
#' @export
#' @seealso \code{\link{scan1perm}}
chr_lengths <-
    function(map, collapse_to_AX=FALSE)
{

    result <- vapply(map, function(a) diff(range(a, na.rm=TRUE)), 1)
    attr(result, "is_x_chr") <- attr(map, "is_x_chr")

    if(collapse_to_AX) return(collapse_chr_lengths_to_AX(result))

    result
}

# turn lengths into sums of autosomal and X-chr lengths
collapse_chr_lengths_to_AX <-
function(lengths, is_x_chr=NULL)
{
    if(length(lengths)==2 && (is.null(names(lengths)) || all(names(lengths)==c("A","X")))) {
        # here, assume that the lengths have already been collapsed
        result <- lengths
        names(result) <- c("A", "X")
        attr(result, "is_x_chr") <- c(A=FALSE,X=TRUE)
        return(result)
    }

    if(is.null(is_x_chr)) # not provided as argument; use attribute
        is_x_chr <- attr(lengths, "is_x_chr")

    if(is.null(is_x_chr)) { # still not provided; assume all autosomes
        warning("No is_x_chr attribute found; assuming all autosomes.")
        is_x_chr <- rep(FALSE, length(lengths))
        names(is_x_chr) <- names(lengths)
    }

    if(length(is_x_chr) != length(lengths)) {
        if(!all(names(lengths) %in% names(is_x_chr)))
           stop("Mismatch between lengths and its is_x_chr attribute")
        is_x_chr <- is_x_chr[names(lengths)]
    }

    result <- c(A=sum(lengths[!is_x_chr]),
                X=sum(lengths[is_x_chr]))
    attr(result, "is_x_chr") <- c(A=FALSE, X=TRUE)

    result
}

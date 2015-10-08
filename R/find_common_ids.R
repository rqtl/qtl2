# match ids
#
# takes two vectors of character strings
# returns a list with
#   first,second = indexes in both vectors, to make them aligned
#   firstonly = indexes in first vector that are not in the second
#   secondonly = indexes in second vector that are not in the first
match_ids <-
    function(x,y)
{
    mxy <- match(x, y)
    myx <- match(y, x)

    result <- list(first=which(!is.na(mxy)),
                   second=mxy[!is.na(mxy)],
                   firstonly=which(is.na(mxy)),
                   secondonly=which(is.na(myx)))

    if(length(result$first) > 0)
        names(result$first) <- x[result$first]
    if(length(result$firstonly) > 0)
        names(result$firstonly) <- x[result$firstonly]
    if(length(result$second) > 0)
        names(result$second) <- y[result$second]
    if(length(result$secondonly) > 0)
        names(result$secondonly) <- y[result$secondonly]

    result
}

# find common ids across a set of vectors of character string
# return value: ids that are present in *all* of the inputs
find_common_ids <-
    function(...)
{
    inp <- list(...)
    n_inp <- length(inp)

    if(n_inp==0)
        stop("No input")

    if(n_inp==1) {
        if(length(inp[[1]])==0) return(character(0))
        return(inp[[1]])
    }

    reduced <- inp[[1]]
    for(i in 2:n_inp) {
        m <- match_ids(reduced, inp[[i]])
        reduced <- reduced[m$first]
    }

    if(length(reduced)==0) return(character(0))

    reduced
}

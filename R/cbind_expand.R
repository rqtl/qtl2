# cbind, matching rownames and expanding with NAs as needed
#
# does cbind(mat1, mat2)
# but: - looks at rownames and makes sure they line up
#      - if rows are in one matrix but not the other,
#        creates row of NAs where missing
cbind_expand <-
    function(...)
{
    input <- list(...)
    if(length(input)<2)
        stop("Must have >= 2 input matrices")

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

# cbind, but allowing different numbers of rows and expanding with NAs as needed
#
# rownames are stripped off
#
# Note that there's a similar function in qtl2geno, but this one is
# different as it doesn't attempt to align the rows by their row names
cbind_expand <-
    function(...)
{
    input <- list(...)
    if(length(input)<2) return(input[[1]])

    if(!all(vapply(input, is.matrix, TRUE)))
        stop("Not all inputs are matrices")

    nrow <- vapply(input, nrow, 1)
    max_nrow <- max(nrow)

    for(i in seq(along=input)) {
        rownames(input[[i]]) <- NULL # strip off the rownames
        if(nrow[i] < max_nrow) { # need to pad
            new_rows <- matrix(NA, ncol=ncol(input[[i]]), nrow=max_nrow - nrow[i])
            colnames(new_rows) <- colnames(input[[i]])
            input[[i]] <- rbind(input[[i]], new_rows)
        }
    }

    do.call("cbind", input)
}

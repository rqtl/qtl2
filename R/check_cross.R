# functions to check the integrity of QTL cross data

# check if a crosstype is supported
#   if give_error=TRUE, problem results in error
#   if give_error=FALSE, returns TRUE if supported; FALSE if not
check_crosstype <-
function(crosstype, give_error=TRUE)
{
    z <- try(.crosstype_supported(crosstype), silent=TRUE)

    if("try-error" %in% class(z) || !z) {
        if(give_error)
            stop("Cross type ", crosstype, " not yet supported.")
        return(invisible(FALSE))
    }

    invisible(TRUE) # if no error, return TRUE
}


# count the number of invalid genotypes
# returns a matrix individuals x chromosomes
count_invalid_genotypes <-
function(cross)
{
    result <- matrix(nrow=nrow(cross$geno[[1]]),
                     ncol=length(cross$geno))
    dimnames(result) <- list(rownames(cross$geno[[1]]),
                             names(cross$geno))

    cross_info <- t(cross$cross_info)

    for(i in seq(along=cross$geno))
        result[,i] <- .count_invalid_genotypes(cross$crosstype,
                                               t(cross$geno[[i]]),
                                               cross$is_x_chr[i],
                                               cross$is_female,
                                               cross_info)

    result
}


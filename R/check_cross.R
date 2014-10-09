# functions to check the integrity of QTL cross data

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


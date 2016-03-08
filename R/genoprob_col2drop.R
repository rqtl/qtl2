# find columns in a genotype probability array that should be dropped
#
# Xonly: If TRUE and probs is a list, just look at the X chromosome
genoprobs_col2drop <-
    function(probs, Xonly=TRUE, tol=1e-8)
{
    if(is.list(probs)) { # proper calc_genoprob object, hopefully
        if(Xonly) {
            result <- lapply(seq(along=probs$probs), function(a) numeric(0))
            names(result) <- probs$chrID

            is_x_chr <- probs$is_x_chr
            if(is.null(is_x_chr) || !any(is_x_chr)) return(result) # no X chr

            w <- which(is_x_chr)
            for(i in w)
                result[[i]] <- genoprobs_col2drop(probs$probs[[i]])
            return(result)
        }
        else {
            return(lapply(probs$probs, genoprobs_col2drop))
        }
    }

    # put genotypes first and take rowMeans
    wh <- which(rowMeans(aperm(probs, c(2,1,3))) < tol)

    names(wh) <- NULL # eliminate names
    wh
}

# find columns in a genotype probability array that should be dropped
genoprobs_col2drop <-
    function(probs, tol=1e-8)
{
    if(is.list(probs))
        return(lapply(probs, genoprobs_col2drop))

    d <- dim(probs)
    # put genotypes first and take rowMeans
    which(rowMeans(aperm(probs, c(2,1,3))) < tol)
}

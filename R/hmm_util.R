# various internal utility functions

# create empty set of matrices for founder genotype data
create_empty_founder_geno <-
function(geno)
{
    result <- vector("list", length(geno))
    names(result) <- names(geno)
    for(i in seq(along=geno)) {
        result[[i]] <- matrix(0L, nrow=0, ncol=ncol(geno[[i]]))
        colnames(result[[i]]) <- colnames(geno[[i]])
    }
    result
}

# convert map to rf
map2rf <-
    function(map, map_function=c("haldane", "kosambi", "c-f", "morgan"), tol=1e-8)
{
    map_function <- match.arg(map_function)

    # if list, use lapply and call this same function for each chromosome
    if(is.list(map))
        return(lapply(map, map2rf, map_function, tol))

    rf <- mf(diff(map), map_function)
    rf[rf < tol] <- tol # don't let values be too small

    rf
}

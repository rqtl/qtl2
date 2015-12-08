# code related to clusters

# test if input is a prepared cluster
is_cluster <-
    function(cores)
{
    "cluster" %in% class(cores) && "SOCKcluster" %in% class(cores)
}

# number of cores being used
n_cores <-
    function(cores)
{
    if(is_cluster(cores)) return( length(cores) )
    cores
}

# set up a cluster
setup_cluster <-
    function(cores, quiet=TRUE)
{
    if(is_cluster(cores)) return(cores)

    if(cores==0) cores <- parallel::detectCores() # if 0, detect cores

    if(cores > 1 && Sys.info()[1] == "Windows") { # windows doesn't support mclapply
        cores <- parallel::makeCluster(cores)
        do.call("on.exit",
                list(quote(parallel::stopCluster(cores))),
                envir=parent.frame())
    }
    cores
}

# run code by cluster
# (to deal with different methods on different architectures)
run_by_cluster <-
    function(cores, ...)
{
    if(is_cluster(cores)) { # cluster object; use mclapply
        return( parallel::clusterApply(cores, ...) )
    } else {
        if(cores==1) return( lapply(...) )
        return( parallel::mclapply(..., mc.cores=cores) )
    }
}

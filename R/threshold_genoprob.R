threshold_genoprob <-
    function(probs, threshold=1e-6, cores=1)
{
    attrib <- attributes(probs)

    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, seq_along(probs),
                             function(i) {
        this_result <- aperm( .threshold_genoprob(aperm(probs[[i]], c(2,1,3)), threshold), c(2,1,3))
        dimnames(this_result) <- dimnames(probs[[i]])
        this_result })

    for(a in names(attrib)) {
        attr(result, a) <- attrib[[a]]
    }

    result
}

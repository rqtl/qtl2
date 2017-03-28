# sim_geno2
# Simulate genotypes given observed marker data
# this version pre-calculates init, step, and emit (attempting to be faster for DO)
#
# Same input and output as sim_geno()
sim_geno2 <-
    function(cross, map=NULL, n_draws=1, error_prob=1e-4,
             map_function=c("haldane", "kosambi", "c-f", "morgan"),
             quiet=TRUE, cores=1)
{
    # check inputs
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)

    # set up cluster; make quiet=FALSE if cores>1
    cores <- setup_cluster(cores)
    if(!quiet && n_cores(cores) > 1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # no more messages
    }

    # pseudomarker map
    if(is.null(map))
        map <- insert_pseudomarkers(cross$gmap)
    # possibly subset the map
    if(length(map) != length(cross$geno) || !all(names(map) == names(cross$geno))) {
        chr <- names(cross$geno)
        if(!all(chr %in% names(map)))
            stop("map doesn't contain all of the necessary chromosomes")
        map <- map[chr]
    }
    # calculate marker index object
    index <- create_marker_index(lapply(cross$geno, colnames), map)

    probs <- vector("list", length(map))
    rf <- map2rf(map, map_function)

    # deal with missing information
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_group_func <- function(i) {
        dr <- .sim_geno2(cross$crosstype, t(cross$geno[[chr]][group[[i]],,drop=FALSE]),
                         founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[group[[i]][1]],
                         cross$cross_info[group[[i]][1],], rf[[chr]], index[[chr]],
                         error_prob, n_draws)
        aperm(dr, c(3,1,2))
    }

    # split individuals into groups with common sex and cross_info
    sex_crossinfo <- paste(cross$is_female, apply(cross$cross_info, 1, paste, collapse=":"), sep=":")
    group <- split(seq(along=sex_crossinfo), sex_crossinfo)
    names(group) <- NULL
    nc <- n_cores(cores)
    while(nc > length(group) && max(sapply(group, length)) > 1) { # successively split biggest group in half until there are as many groups as cores
        mx <- which.max(sapply(group, length))
        g <- group[[mx]]
        group <- c(group, list(g[seq(1, length(g), by=2)]))
        group[[mx]] <- g[seq(2, length(g), by=2)]
    }
    groupindex <- seq(along=group)

    draws <- vector("list", length(cross$geno))
    names(draws) <- names(cross$geno)
    for(chr in seq(along=cross$geno)) {
        if(!quiet) message("Chr ", names(cross$geno)[chr])

        # calculations in parallel [if cores==1, it just does lapply()]
        temp <- cluster_lapply(cores, groupindex, by_group_func)

        # paste them back together
        d <- vapply(temp, dim, rep(0,3))
        nr <- sum(d[1,])
        draws[[chr]] <- array(dim=c(nr, d[2,1], d[3,1]))
        for(i in groupindex)
            draws[[chr]][group[[i]],,] <- temp[[i]]

        dimnames(draws[[chr]]) <- list(rownames(cross$geno[[chr]]),
                                       names(map[[chr]]),
                                       NULL)

    }

    names(draws) <- names(cross$gmap)

    attr(draws, "crosstype") <- cross$crosstype
    attr(draws, "is_x_chr") <- cross$is_x_chr
    attr(draws, "alleles") <- cross$alleles

    class(draws) <- c("sim_geno", "list")
    draws
}

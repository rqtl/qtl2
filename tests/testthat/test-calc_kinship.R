context("Calculation of kinship matrix")

test_that("calc_kinship works for RIL", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    grid <- calc_grid(grav2$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)
    sim <- calc_kinship(probs_sub)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(grav2$geno[[1]]))
    expect_equal(colnames(sim), rownames(grav2$geno[[1]]))

    # check a few pairs
    set.seed(88213118)
    pairs <- list(c(1,2))
    n <- n_ind(grav2)
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n, 2), sample(n, 2)))
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    for(i in seq(along=probs)) {
        for(k in seq(along=pairs))
            expected[k] <- expected[k] + sum(probs_sub[[i]][pairs[[k]][1],,] * probs_sub[[i]][pairs[[k]][2],,])
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs))
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])

})

test_that("calc_kinship works for F2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    grid <- calc_grid(iron$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)
    sim <- calc_kinship(probs_sub, omit_x=TRUE)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(iron$geno[[1]]))
    expect_equal(colnames(sim), rownames(iron$geno[[1]]))

    f2_geno2alle <-
        function(prob, x_chr=FALSE)
        {
            if(x_chr) {
                prob[,1,] <- prob[,1,]+prob[,2,]/2+prob[,3,]/2+prob[,5,]
                prob[,2,] <- prob[,4,]+prob[,2,]/2+prob[,3,]/2+prob[,6,]
                return(prob[,1:2,])

            } else {
                prob[,1,] <- prob[,1,]+prob[,2,]/2
                prob[,2,] <- prob[,3,]+prob[,2,]/2
                return(prob[,1:2,])
            }
        }

    # check a few values
    set.seed(54028069)
    n <- nrow(probs[[1]])
    pairs <- list(c(1,2))
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n, 2), sample(n, 2)))
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    is_x_chr <- attr(probs_sub, "is_x_chr")
    for(i in which(!is_x_chr)) { # just use autosomes
        for(k in seq(along=pairs)) {
            prob_1 <- probs_sub[[i]][pairs[[k]][1],,,drop=FALSE]
            prob_2 <- probs_sub[[i]][pairs[[k]][2],,,drop=FALSE]
            # autosomes: convert to allele probs
            prob_1 <- f2_geno2alle(prob_1)
            prob_2 <- f2_geno2alle(prob_2)

            expected[k] <- expected[k] + sum(prob_1 * prob_2)
        }
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs))
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])

    # version using genotype probabilities
    sim <- calc_kinship(probs_sub, use_allele_probs=FALSE, omit_x=TRUE)
    sim2 <- calc_kinship(probs_sub, use_allele_probs=FALSE, omit_x=TRUE)
    expect_equal(sim, sim2)

    # check a few values
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    for(i in which(!is_x_chr)) { # just use autosomes
        for(k in seq(along=pairs)) {
            prob_1 <- probs_sub[[i]][pairs[[k]][1],,,drop=FALSE]
            prob_2 <- probs_sub[[i]][pairs[[k]][2],,,drop=FALSE]

            expected[k] <- expected[k] + sum(prob_1 * prob_2)
        }
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs))
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])


    # also try with X chr (now the default)
    sim <- calc_kinship(probs_sub)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(iron$geno[[1]]))
    expect_equal(colnames(sim), rownames(iron$geno[[1]]))

    # check a few values
    set.seed(54028069)
    n <- nrow(probs[[1]])
    pairs <- list(c(1,2))
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n, 2), sample(n, 2)))
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    for(i in seq(along=probs_sub)) {
        for(k in seq(along=pairs)) {
            prob_1 <- probs_sub[[i]][pairs[[k]][1],,,drop=FALSE]
            prob_2 <- probs_sub[[i]][pairs[[k]][2],,,drop=FALSE]
            prob_1 <- f2_geno2alle(prob_1, is_x_chr[i])
            prob_2 <- f2_geno2alle(prob_2, is_x_chr[i])

            expected[k] <- expected[k] + sum(prob_1 * prob_2)
        }
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs)) {
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])
    }

    # version using genotype probabilities
    sim <- calc_kinship(probs_sub, omit_x=FALSE, use_allele_probs=FALSE)
    sim2 <- calc_kinship(probs_sub, omit_x=FALSE, use_allele_probs=FALSE)
    expect_equal(sim, sim2)

    # check a few values
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    for(i in seq(along=probs_sub)) {
        for(k in seq(along=pairs)) {
            prob_1 <- probs_sub[[i]][pairs[[k]][1],,,drop=FALSE]
            prob_2 <- probs_sub[[i]][pairs[[k]][2],,,drop=FALSE]

            expected[k] <- expected[k] + sum(prob_1 * prob_2)
        }
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs))
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])

})

test_that("calc_kinship chr & loco work for F2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    grid <- calc_grid(iron$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)
    sim <- calc_kinship(probs_sub, omit_x=TRUE)

    sim_chr <- calc_kinship(probs_sub, "chr", omit_x=TRUE)
    sim_loco <- calc_kinship(probs_sub, "loco", omit_x=TRUE)

    # combine results from sim_chr and compare to sim
    sim_combchr <- sim_chr[[1]]*attr(sim_chr[[1]], "n_pos")
    for(i in 2:19)
        sim_combchr <- sim_combchr + sim_chr[[i]]*attr(sim_chr[[i]], "n_pos")
    totpos <- sum(sapply(sim_chr[1:19], attr, "n_pos"))
    sim_combchr <- sim_combchr/totpos
    attr(sim_combchr, "n_pos") <- totpos
    expect_equal(sim_combchr, sim)

    # calculate results with one chromosome at a time
    sim_alt <- vector("list", length(probs))
    names(sim_alt) <- names(probs)
    for(i in names(probs))
        sim_alt[[i]] <- calc_kinship(probs_sub[,i], omit_x=FALSE)
    expect_equal(sim_alt, sim_chr)

    # compare sim - sim_chr with sim_loco
    sim_loco_alt <- vector("list", length(probs))
    names(sim_loco_alt) <- names(probs)
    is_x_chr <- attr(probs, "is_x_chr")
    totpos <- attr(sim, "n_pos")
    for(i in seq(along=probs)) {
        if(is_x_chr[i])
            sim_loco_alt[[i]] <- sim
        else {
            npos <- attr(sim_chr[[i]], "n_pos")
            sim_loco_alt[[i]] <- (sim*totpos - sim_chr[[i]]*npos)/(totpos-npos)
            attr(sim_loco_alt[[i]], "n_pos") <- totpos-npos
        }
    }
    expect_equal(sim_loco_alt, sim_loco)

})

test_that("calc_kinship chr & loco work for F2, when including X", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    grid <- calc_grid(iron$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)
    sim <- calc_kinship(probs_sub, omit_x=FALSE) # *not* omitting X chr

    sim_chr <- calc_kinship(probs_sub, "chr")
    sim_loco <- calc_kinship(probs_sub, "loco", omit_x=FALSE) # *not* omitting X chr

    # combine results from sim_chr and compare to sim
    sim_combchr <- sim_chr[[1]]*attr(sim_chr[[1]], "n_pos")
    for(i in 2:20)
        sim_combchr <- sim_combchr + sim_chr[[i]]*attr(sim_chr[[i]], "n_pos")
    totpos <- sum(sapply(sim_chr, attr, "n_pos"))
    sim_combchr <- sim_combchr/totpos
    attr(sim_combchr, "n_pos") <- totpos
    expect_equal(sim_combchr, sim)

    # compare sim - sim_chr with sim_loco
    sim_loco_alt <- vector("list", length(probs))
    names(sim_loco_alt) <- names(probs)
    totpos <- attr(sim, "n_pos")
    for(i in seq(along=probs)) {
        npos <- attr(sim_chr[[i]], "n_pos")
        sim_loco_alt[[i]] <- (sim*totpos - sim_chr[[i]]*npos)/(totpos-npos)
        attr(sim_loco_alt[[i]], "n_pos") <- totpos-npos
    }
    expect_equal(sim_loco_alt, sim_loco)

})


test_that("calc_kinship chr & loco work when multi-core", {
    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    grid <- calc_grid(iron$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)

    sim_chr <- calc_kinship(probs_sub, "chr")
    sim_chr_mc <- calc_kinship(probs_sub, "chr", cores=4)
    expect_equal(sim_chr_mc, sim_chr)

    sim_loco <- calc_kinship(probs_sub, "loco")
    sim_loco_mc <- calc_kinship(probs_sub, "loco", cores=4)
    expect_equal(sim_loco_mc, sim_loco)

    sim_loco <- calc_kinship(probs_sub, "loco", omit_x=FALSE)
    sim_loco_mc <- calc_kinship(probs_sub, "loco", omit_x=FALSE, cores=4)
    expect_equal(sim_loco_mc, sim_loco)


    ## RIL
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    grid <- calc_grid(grav2$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)

    sim_chr <- calc_kinship(probs_sub, "chr")
    sim_chr_mc <- calc_kinship(probs_sub, "chr", cores=4)
    expect_equal(sim_chr_mc, sim_chr)

    sim_loco <- calc_kinship(probs_sub, "loco")
    sim_loco_mc <- calc_kinship(probs_sub, "loco", cores=4)
    expect_equal(sim_loco_mc, sim_loco)

})

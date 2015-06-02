context("Calculation of genetic similarity")

test_that("calc_genetic_sim works for RIL", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    sim <- calc_genetic_sim(probs)

    # pre-subset to grid
    probs_sub <- probs_to_grid(probs)
    sim2 <- calc_genetic_sim(probs_sub)
    expect_equal(sim, sim2)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(grav2$geno[[1]]))
    expect_equal(colnames(sim), rownames(grav2$geno[[1]]))

    # check a few values
    set.seed(88213118)
    n_ind <- nrow(probs[[1]])
    pairs <- list(c(1,2))
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n_ind, 2), sample(n_ind, 2)))
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

test_that("calc_genetic_sim works for F2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)
    sim <- calc_genetic_sim(probs)

    # pre-subset to grid
    probs_sub <- probs_to_grid(probs)
    sim2 <- calc_genetic_sim(probs_sub)
    expect_equal(sim, sim2)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(iron$geno[[1]]))
    expect_equal(colnames(sim), rownames(iron$geno[[1]]))

    f2_geno2alle <-
        function(prob)
        {
            prob[,1,] <- prob[,1,]+prob[,2,]/2
            prob[,2,] <- prob[,3,]+prob[,2,]/2
            prob[,1:2,]
        }

    # check a few values
    set.seed(54028069)
    n_ind <- nrow(probs[[1]])
    pairs <- list(c(1,2))
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n_ind, 2), sample(n_ind, 2)))
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
    sim <- calc_genetic_sim(probs, use_allele_probs=FALSE)
    sim2 <- calc_genetic_sim(probs, use_allele_probs=FALSE)
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


    # also try with X chr
    sim <- calc_genetic_sim(probs, omit_x=FALSE)

    # pre-subset to grid
    sim2 <- calc_genetic_sim(probs_sub, omit_x=FALSE)
    expect_equal(sim, sim2)

    # row and colnames okay
    expect_equal(rownames(sim), rownames(iron$geno[[1]]))
    expect_equal(colnames(sim), rownames(iron$geno[[1]]))

    # check a few values
    set.seed(54028069)
    n_ind <- nrow(probs[[1]])
    pairs <- list(c(1,2))
    for(i in 1:5)
        pairs <- c(pairs, list(sample(n_ind, 2), sample(n_ind, 2)))
    expected <- rep(0, length(pairs))
    tot_pos <- 0
    for(i in seq(along=probs_sub)) {
        for(k in seq(along=pairs)) {
            prob_1 <- probs_sub[[i]][pairs[[k]][1],,,drop=FALSE]
            prob_2 <- probs_sub[[i]][pairs[[k]][2],,,drop=FALSE]
            if(!is_x_chr[i]) { # autosomes: convert to allele probs
                prob_1 <- f2_geno2alle(prob_1)
                prob_2 <- f2_geno2alle(prob_2)
            }

            expected[k] <- expected[k] + sum(prob_1 * prob_2)
        }
        tot_pos <- tot_pos + dim(probs_sub[[i]])[3]
    }
    expected <- expected/tot_pos
    for(k in seq(along=pairs))
        expect_equal(sim[pairs[[k]][1],pairs[[k]][2]], expected[k])

    # version using genotype probabilities
    sim <- calc_genetic_sim(probs, omit_x=FALSE, use_allele_probs=FALSE)
    sim2 <- calc_genetic_sim(probs, omit_x=FALSE, use_allele_probs=FALSE)
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

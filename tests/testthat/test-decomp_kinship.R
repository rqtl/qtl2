context("eigen decomposition of kinship matrix")

test_that("eigen decomposition works", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)

    K <- calc_kinship(probs)
    Ke <- decomp_kinship(K)
    expected <- eigen(K)

    expect_equal(sort(Ke$values), sort(expected$values))
    expect_equivalent(K, t(Ke$vectors) %*% diag(Ke$values) %*% Ke$vectors)
    expect_equivalent(solve(K), t(Ke$vectors) %*% diag(1/Ke$values) %*% Ke$vectors)

    # run it through again and get same answer
    Ke2 <- decomp_kinship(Ke)
    expect_equal(Ke2, Ke)

    # list of matrices
    K <- calc_kinship(probs, "loco")
    Ke <- decomp_kinship(K)

    for(i in seq(along=K)) {
        expected <- eigen(K[[i]])

        expect_equal(sort(Ke[[i]]$values), sort(expected$values))
        expect_equivalent(K[[i]], t(Ke[[i]]$vectors) %*% diag(Ke[[i]]$values) %*% Ke[[i]]$vectors)
        expect_equivalent(solve(K[[i]]), t(Ke[[i]]$vectors) %*% diag(1/Ke[[i]]$values) %*% Ke[[i]]$vectors)
    }

    # run it through again and get same answer
    Ke2 <- decomp_kinship(Ke)
    expect_equal(Ke2, Ke)

})

test_that("multi-core eigen decomposition works", {
    if(isnt_karl()) skip("this test only run locally")

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)

    K <- calc_kinship(probs, "loco", cores=4)
    Ke <- decomp_kinship(K)
    Ke_multicore <- decomp_kinship(K, cores=4)
    expect_equal(Ke_multicore, Ke)

    K <- calc_kinship(probs, "chr", cores=4)
    Ke <- decomp_kinship(K)
    Ke_multicore <- decomp_kinship(K, cores=4)
    expect_equal(Ke_multicore, Ke)

})

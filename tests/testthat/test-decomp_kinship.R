context("eigen decomposition of kinship matrix")

test_that("eigen decomposition works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron[,c(1,3,5)], map, error_prob=0.002) # <- use just 3 chromosomes, for speed

    K <- calc_kinship(probs)
    Ke <- decomp_kinship(K)
    expected <- eigen(K)

    expect_equal(sort(Ke$values), sort(expected$values))
    expect_equivalent(K, t(Ke$vectors) %*% diag(Ke$values) %*% Ke$vectors)

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
    }

    # run it through again and get same answer
    Ke2 <- decomp_kinship(Ke)
    expect_equal(Ke2, Ke)

    # make sure it gives an error input has 0 rows and 0 cols
    expect_error(decomp_kinship(K[[1]][numeric(0),numeric(0),drop=FALSE]))

    # make sure it gives an error input isn't square
    expect_error(decomp_kinship(K[[1]][1:5,1:4]))

    # make sure it gives an error input isn't a matrix
    expect_error(decomp_kinship(K[[1]][5,5]))

})

test_that("multi-core eigen decomposition works", {
    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    K <- calc_kinship(probs, "loco", cores=4)
    Ke <- decomp_kinship(K)
    Ke_multicore <- decomp_kinship(K, cores=4)
    expect_equal(Ke_multicore, Ke)

    K <- calc_kinship(probs, "chr", cores=4)
    Ke <- decomp_kinship(K)
    Ke_multicore <- decomp_kinship(K, cores=4)
    expect_equal(Ke_multicore, Ke)

})


test_that("eigen decomposition gives error with non-square matrix", {

    k <- matrix(runif(100), ncol=5)
    expect_error(decomp_eigen(k))

    k <- list(k, matrix(runif(100), ncol=20))
    expect_error(decomp_eigen(k))

})

test_that("eigen decomposition works with scan1", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron[,c(1,3,5)], map, error_prob=0.002) # <- use just 3 chromosomes, for speed

    K <- calc_kinship(probs)
    Ke <- decomp_kinship(K)

    out <- scan1(probs, iron$pheno, K)
    expect_equal(scan1(probs, iron$pheno, Ke), out)

    K_loco <- calc_kinship(probs, "loco")
    Ke_loco <- decomp_kinship(K_loco)

    out_loco <- scan1(probs, iron$pheno, K_loco)
    expect_equal(scan1(probs, iron$pheno, Ke_loco), out_loco)

})

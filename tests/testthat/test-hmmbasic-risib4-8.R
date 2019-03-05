context("basic HMM functions in 4- and 8-way RIL by sib-mating")

test_that("risib4-8 n_gen, n_alleles work", {

    expect_equal(nalleles("risib4"), 4)
    expect_equal(nalleles("risib8"), 8)

    expect_equal(test_ngen("risib4", FALSE), 4)
    expect_equal(test_ngen("risib8", FALSE), 8)

    # X chromosome
    expect_equal(test_ngen("risib4", TRUE), 4)
    expect_equal(test_ngen("risib8", TRUE), 8)

})

test_that("risib4-8 possible_gen work", {

    expect_equal(test_possible_gen("risib4", FALSE, FALSE, 1:4), 1:4)
    expect_equal(test_possible_gen("risib8", FALSE, FALSE, 1:8), 1:8)

    # X chr
    expect_equal(test_possible_gen("risib4", TRUE, FALSE, 1:4), 1:3)
    expect_equal(test_possible_gen("risib8", TRUE, FALSE, 1:8), c(1:3,5:6))

    # X chr, different cross order
    expect_equal(test_possible_gen("risib4", TRUE, FALSE, c(3,4,1,2)), c(3,4,1))
    expect_equal(test_possible_gen("risib8", TRUE, FALSE, c(4,5,1,3,2,7,8,6)), c(4,5,1,2,7))

    # one more time
    expect_equal(test_possible_gen("risib4", TRUE, FALSE, c(2,1,4,3)), c(2,1,4))
    expect_equal(test_possible_gen("risib8", TRUE, FALSE, c(1,5,6,8,4,2,7,3)), c(1,5,6,4,2))

})

# FIX_ME: test check_geno

test_that("risib4-8 init work", {

    expect_equal( sapply(1:4, function(i) test_init("risib4", i, FALSE, FALSE, 1:4)), log(rep(1/4, 4)) )
    expect_equal( sapply(1:8, function(i) test_init("risib8", i, FALSE, FALSE, 1:8)), log(rep(1/8, 8)) )

    # X chromosome
    expect_equal( sapply(1:3, function(i) test_init("risib4", i, TRUE, FALSE, 1:4)), log(rep(1/3, 3)) )
    expect_equal( sapply(c(1:3,5,6), function(i) test_init("risib8", i, TRUE, FALSE, 1:8)), log(c(1,1,2,1,1)/6) )

    # X chromosoem, different order
    expect_equal( sapply(c(3,4,1), function(i) test_init("risib4", i, TRUE, FALSE, c(3,4,1,2))), log(rep(1/3, 3)) )
    expect_equal( sapply(c(4,5,1,2,7), function(i) test_init("risib8", i, TRUE, FALSE, c(4,5,1,3,2,7,8,6))), log(c(1,1,2,1,1)/6) )

    # one more time
    expect_equal( sapply(c(1,2,4), function(i) test_init("risib4", i, TRUE, FALSE, c(2,1,4,3))), log(rep(1/3, 3)) )
    expect_equal( sapply(c(1,5,6,4,2), function(i) test_init("risib8", i, TRUE, FALSE, c(1,5,6,8,4,2,7,3))), log(c(1,1,2,1,1)/6) )

})

# FIX_ME: test emit

test_that("risib4 step works for autosome", {

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf*2/(1+6*rf), ncol=4, nrow=4)
        diag(expected) <- 1/(1+6*rf)

        result <- matrix(ncol=4, nrow=4)
        for(i in 1:4) {
            for(j in 1:4) {
                result[i,j] <- test_step("risib4", i, j, rf, FALSE, FALSE, 1:4)
            }
        }

        expect_equal(result, log(expected))
    }

    # test_stepmatrix: output has rows that sum to 1
    obs <- vapply(test_stepmatrix("risib4", c(0.01, 0.1, 0.45), FALSE, FALSE, 1:4),
                  function(a) rowSums(exp(a)), rep(1,4))
    expect_equal(obs, matrix(1, ncol=3, nrow=4))

    obs <- vapply(test_stepmatrix("risib4", c(0.01, 0.1, 0.45), FALSE, FALSE, c(3,1,2,4)),
                  function(a) rowSums(exp(a)), rep(1,4))
    expect_equal(obs, matrix(1, ncol=3, nrow=4))

    obs <- vapply(test_stepmatrix("risib4", c(0.01, 0.1, 0.45), TRUE, FALSE, c(3,1,2,4)),
                  function(a) rowSums(exp(a)), rep(1,3))
    expect_equal(obs, matrix(1, ncol=3, nrow=3))

})


test_that("risib8 step works for autosome", {

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/(1+6*rf), ncol=8, nrow=8)
        diag(expected) <- (1-rf)/(1+6*rf)

        result <- matrix(ncol=8, nrow=8)
        for(i in 1:8) {
            for(j in 1:8) {
                result[i,j] <- test_step("risib8", i, j, rf, FALSE, FALSE, 1:8)
            }
        }

        expect_equal(result, log(expected))
    }

    # test_stepmatrix: output has rows that sum to 1
    obs <- vapply(test_stepmatrix("risib8", c(0.01, 0.1, 0.45), FALSE, FALSE, 1:8),
                  function(a) rowSums(exp(a)), rep(1,8))
    expect_equal(obs, matrix(1, ncol=3, nrow=8))

    obs <- vapply(test_stepmatrix("risib8", c(0.01, 0.1, 0.45), FALSE, FALSE, c(3,1,5:8,2,4)),
                  function(a) rowSums(exp(a)), rep(1,8))
    expect_equal(obs, matrix(1, ncol=3, nrow=8))

    obs <- vapply(test_stepmatrix("risib8", c(0.01, 0.1, 0.45), TRUE, FALSE, c(3,1,5:8,2,4)),
                  function(a) rowSums(exp(a)), rep(1,5))
    expect_equal(obs, matrix(1, ncol=3, nrow=5))

})



test_that("risib4-8 geno_names work", {

    expect_equal( geno_names("risib4", LETTERS[5:8], FALSE), paste0(LETTERS[5:8], LETTERS[5:8]) )
    expect_equal( geno_names("risib8", LETTERS[2:9], FALSE), paste0(LETTERS[2:9], LETTERS[2:9]) )

})

test_that("risib4-8 nrec work", {

    x <- matrix(ncol=8, nrow=8)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=8)

    res4 <- matrix(ncol=4, nrow=4)
    res8 <- matrix(ncol=8, nrow=8)
    for(i in 1:8) {
        for(j in 1:8) {
            if(i<5 && j<5) res4[i,j] <- test_nrec("risib4", i, j, FALSE, FALSE, 1:4)
            if(i<9 && j<9) res8[i,j] <- test_nrec("risib8", i, j, FALSE, FALSE, 1:8)
        }
    }

    expect_equal( res4, x[1:4,1:4])
    expect_equal( res8, x)

})

context("basic HMM functions in Heterogeneous Stock F1 (crossed to a 9th inbred strain)")

test_that("HSF1 nalleles works", {
    expect_equal(nalleles("hsf1"), 8)
})

test_that("HSF1 check_geno works", {

    # observed genotypes
    for(i in 0:5) {
        # Autosome
        expect_true(test_check_geno("hsf1", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("hsf1", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_true(test_check_geno("hsf1", i, TRUE, TRUE, FALSE, 20))
    }
    for(i in c(-1, 6)) {
        # Autosome
        expect_false(test_check_geno("hsf1", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("hsf1", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_false(test_check_geno("hsf1", i, TRUE, TRUE, FALSE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 1:8) {
        # Autosome
        expect_true(test_check_geno("hsf1", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("hsf1", i, FALSE, TRUE, TRUE, 20))
    }
    for(i in c(0, 9)) {
        # Autosome
        expect_false(test_check_geno("hsf1", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("hsf1", i, FALSE, TRUE, TRUE, 20))
    }

    # true genotypes, male X
    for(i in c(0, 8 + 1:8)) {
        expect_false(test_check_geno("hsf1", i, FALSE, TRUE, FALSE, 20))
    }
    for(i in 1:8) {
        expect_true(test_check_geno("hsf1", i, FALSE, TRUE, FALSE, 20))
    }

})

test_that("HSF1 n_gen works", {

    expect_equal(test_ngen("hsf1", FALSE), 8)
    expect_equal(test_ngen("hsf1", TRUE),  8)

})

test_that("HSF1 possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("hsf1", FALSE, FALSE, 20), 1:8)

    # X female
    expect_equal(test_possible_gen("hsf1", TRUE, TRUE, 20), 1:8)

    # X male
    expect_equal(test_possible_gen("hsf1", TRUE, FALSE, 20), 1:8)

})

test_that("HSF1 init works", {

    hom <- 1:8
    male <- 8 + (1:8)

    for(i in hom) {
        # autosome and female X
        expect_equal(test_init("hsf1", i, FALSE, FALSE, 20), log(1/8))
        expect_equal(test_init("hsf1", i, TRUE,  TRUE,  20), log(1/8))
    }

    for(i in male)
        expect_equal(test_init("hsf1", i, TRUE,  FALSE, 20), log(1/8))

})


test_that("HSF1 emit works", {

    fgen <- c(1,3,0,1,3,0,1,3,1) # founder genotypes; 0=missing, 1=AA, 3=BB
    err <- 0.01

    expected <- log(c(1-err, err/2, err/2, 1-err/2, err))
    for(trueg in c(1,4,7)) { # A, D, G
        for(obsg in 1:5) {
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, TRUE, FALSE, 20), expected[obsg])
        }
    }
    expected <- log(c(err/2, 1-err, err/2, 1-err/2, 1-err/2))
    expected_male <- log(c(err/2, err/2, 1-err, err, 1-err/2))
    for(trueg in c(2,5,8)) { # B, E, H
        for(obsg in 1:5) {
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, TRUE, FALSE, 20), expected_male[obsg])
        }
    }
    expected <- log(c(1-err, 1, err, 1-err, err))
    expected_male <- log(c(1, 1, 1, 1, 1))
    for(trueg in c(3,6)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
            expect_equal(test_emit("hsf1", obsg, trueg, err, fgen, TRUE, FALSE, 20), expected_male[obsg])
        }
    }

})

test_that("HSF1 step works", {

    ng <- 8
    trmat <- matrix(nrow=ng, ncol=ng)
    # autosome
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("hsf1", gl, gr, rf, FALSE, FALSE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))

            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

    # female X
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("hsf1", gl, gr, rf, TRUE, TRUE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))

            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

    # male X
    ng <- 8
    trmat <- matrix(nrow=ng, ncol=ng)
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("hsf1", 8+gl, 8+gr, rf, TRUE, FALSE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))

            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

})

test_that("geno_names works", {
    auto <- LETTERS[1:8]

    expect_equal(geno_names("hsf1", LETTERS[1:8], FALSE), auto)
    expect_equal(geno_names("hsf1", LETTERS[1:8], TRUE), auto)
})

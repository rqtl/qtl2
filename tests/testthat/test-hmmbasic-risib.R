
context("basic HMM functions in risib")

test_that("risib check_geno works", {

    # autosome
    expect_true(test_check_geno("risib", 0, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("risib", 1, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("risib", 2, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("risib", 1, FALSE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("risib", 2, FALSE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("risib", 3, TRUE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("risib", 0, FALSE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("risib", 3, FALSE, FALSE, FALSE, 0, FALSE))

    # X chromosome
    expect_true(test_check_geno("risib", 0, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("risib", 1, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("risib", 2, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("risib", 1, FALSE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("risib", 2, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("risib", 3, TRUE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("risib", 0, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("risib", 3, FALSE, TRUE, TRUE, 0, FALSE))

    # X chr reverse
    expect_true(test_check_geno("risib", 0, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("risib", 1, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("risib", 2, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("risib", 1, FALSE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("risib", 2, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("risib", 3, TRUE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("risib", 0, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("risib", 3, FALSE, TRUE, TRUE, 1, FALSE))

})

test_that("risib all_geno works", {
    expect_equal(test_allgeno("risib", FALSE, FALSE), 1:2)
    expect_equal(test_allgeno("risib", TRUE, FALSE), 1:2)

    # phase-known shouldn't matter
    expect_equal(test_allgeno("risib", FALSE, TRUE), 1:2)
    expect_equal(test_allgeno("risib", TRUE, TRUE), 1:2)
})

test_that("risib geno works", {

    # autosome
    expect_equal(test_geno("risib", FALSE, FALSE, 0, FALSE), 1:2)
    # X chr
    expect_equal(test_geno("risib", TRUE, TRUE, 0, FALSE), 1:2)
    # X chr reverse
    expect_equal(test_geno("risib", TRUE, TRUE, 1, FALSE), 1:2)

    # phase-known shouldn't matter
    expect_equal(test_geno("risib", FALSE, FALSE, 0, TRUE), 1:2)
    expect_equal(test_geno("risib", TRUE, TRUE, 0, TRUE), 1:2)
    expect_equal(test_geno("risib", TRUE, TRUE, 1, TRUE), 1:2)
})

test_that("risib nrec works", {

    # autosome
    expect_equal(test_nrec("risib", 1, 1, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("risib", 1, 2, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("risib", 2, 1, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("risib", 2, 2, FALSE, FALSE, 0), 0)
    # X chr
    expect_equal(test_nrec("risib", 1, 1, TRUE, TRUE, 0), 0)
    expect_equal(test_nrec("risib", 1, 2, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("risib", 2, 1, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("risib", 2, 2, TRUE, TRUE, 0), 0)
    # X chr reverse
    expect_equal(test_nrec("risib", 1, 1, TRUE, TRUE, 1), 0)
    expect_equal(test_nrec("risib", 1, 2, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("risib", 2, 1, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("risib", 2, 2, TRUE, TRUE, 1), 0)

    # some errors
    expect_error(test_nrec("risib", 0, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("risib", 1, 0, FALSE, FALSE, 0))
    expect_error(test_nrec("risib", 3, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("risib", 2, 3, FALSE, FALSE, 0))
    expect_error(test_nrec("risib", 1, 0, TRUE, TRUE, 0))
    expect_error(test_nrec("risib", 0, 2, TRUE, TRUE, 0))
    expect_error(test_nrec("risib", 3, 1, TRUE, TRUE, 0))
    expect_error(test_nrec("risib", 2, 3, TRUE, TRUE, 0))
    expect_error(test_nrec("risib", 1, 0, TRUE, TRUE, 1))
    expect_error(test_nrec("risib", 0, 2, TRUE, TRUE, 1))
    expect_error(test_nrec("risib", 3, 1, TRUE, TRUE, 1))
    expect_error(test_nrec("risib", 2, 3, TRUE, TRUE, 1))

})

test_that("risib init works", {

    # autosome
    expect_equal(test_init("risib", 1, FALSE, FALSE, 0), log(0.5))
    expect_equal(test_init("risib", 2, FALSE, FALSE, 0), log(0.5))
    # X chr
    expect_equal(test_init("risib", 1, TRUE, TRUE, 0), log(2/3))
    expect_equal(test_init("risib", 2, TRUE, TRUE, 0), log(1/3))
    # X chr reverse
    expect_equal(test_init("risib", 1, TRUE, TRUE, 1), log(1/3))
    expect_equal(test_init("risib", 2, TRUE, TRUE, 1), log(2/3))

    # phase-known shouldn't matter
    expect_equal(test_init("risib", 1, FALSE, FALSE, 0, TRUE), log(0.5))
    expect_equal(test_init("risib", 2, FALSE, FALSE, 0, TRUE), log(0.5))
    expect_equal(test_init("risib", 1, TRUE, TRUE, 0, TRUE), log(2/3))
    expect_equal(test_init("risib", 2, TRUE, TRUE, 0, TRUE), log(1/3))
    expect_equal(test_init("risib", 1, TRUE, TRUE, 1, TRUE), log(1/3))
    expect_equal(test_init("risib", 2, TRUE, TRUE, 1, TRUE), log(2/3))

    # errors
    expect_error(test_init("risib", 0, FALSE, FALSE, 0))
    expect_error(test_init("risib", 3, FALSE, FALSE, 0))
    expect_error(test_init("risib", 0, TRUE, TRUE, 0))
    expect_error(test_init("risib", 3, TRUE, TRUE, 0))
    expect_error(test_init("risib", 0, TRUE, FALSE, 1))
    expect_error(test_init("risib", 3, TRUE, FALSE, 1))

    expect_error(test_init("risib", 0, FALSE, FALSE, 0, TRUE))
    expect_error(test_init("risib", 3, FALSE, FALSE, 0, TRUE))
    expect_error(test_init("risib", 0, TRUE, TRUE, 0, TRUE))
    expect_error(test_init("risib", 3, TRUE, TRUE, 0, TRUE))
    expect_error(test_init("risib", 0, TRUE, FALSE, 1, TRUE))
    expect_error(test_init("risib", 3, TRUE, FALSE, 1, TRUE))

})

test_that("risib emit works", {

    # autosome
    eps <- 0.01
    expect_equal(test_emit("risib", 0, 1, eps, FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 0, 2, eps, FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 1, 1, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, FALSE, FALSE, 0), log(1-eps))
    # X chr
    expect_equal(test_emit("risib", 0, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 0, 2, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 1, 1, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, TRUE, TRUE, 0), log(1-eps))
    # X chr reverse
    expect_equal(test_emit("risib", 0, 1, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 0, 2, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 1, 1, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, TRUE, FALSE, 1), log(1-eps))

    # phase-known shouldn't matter
    eps <- 0.001
    expect_equal(test_emit("risib", 0, 1, eps, FALSE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("risib", 0, 2, eps, FALSE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("risib", 1, 1, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, FALSE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, FALSE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    # X chr
    expect_equal(test_emit("risib", 0, 1, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("risib", 0, 2, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("risib", 1, 1, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    # X chr reverse
    expect_equal(test_emit("risib", 0, 1, eps, TRUE, FALSE, 1, TRUE), 0)
    expect_equal(test_emit("risib", 0, 2, eps, TRUE, FALSE, 1, TRUE), 0)
    expect_equal(test_emit("risib", 1, 1, eps, TRUE, FALSE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, TRUE, FALSE, 1, TRUE), log(1-eps))

    # errors
    expect_error(test_emit("risib", 0, 0, eps, FALSE, FALSE, 0))
    expect_error(test_emit("risib", 0, 3, eps, FALSE, FALSE, 0))
    expect_error(test_emit("risib", 3, 1, eps, FALSE, FALSE, 0))
    # X chr
    expect_error(test_emit("risib", 0, 0, eps, TRUE, TRUE, 0))
    expect_error(test_emit("risib", 0, 3, eps, TRUE, TRUE, 0))
    expect_error(test_emit("risib", 3, 1, eps, TRUE, TRUE, 0))
    # X chr reverse
    expect_error(test_emit("risib", 0, 0, eps, TRUE, FALSE, 1))
    expect_error(test_emit("risib", 0, 3, eps, TRUE, FALSE, 1))
    expect_error(test_emit("risib", 3, 1, eps, TRUE, FALSE, 1))

    # errors
    expect_error(test_emit("risib", 0, 0, eps, FALSE, FALSE, 0, TRUE))
    expect_error(test_emit("risib", 0, 3, eps, FALSE, FALSE, 0, TRUE))
    expect_error(test_emit("risib", 3, 1, eps, FALSE, FALSE, 0, TRUE))
    # X chr
    expect_error(test_emit("risib", 0, 0, eps, TRUE, TRUE, 0, TRUE))
    expect_error(test_emit("risib", 0, 3, eps, TRUE, TRUE, 0, TRUE))
    expect_error(test_emit("risib", 3, 1, eps, TRUE, TRUE, 0, TRUE))
    # X chr reverse
    expect_error(test_emit("risib", 0, 0, eps, TRUE, FALSE, 1, TRUE))
    expect_error(test_emit("risib", 0, 3, eps, TRUE, FALSE, 1, TRUE))
    expect_error(test_emit("risib", 3, 1, eps, TRUE, FALSE, 1, TRUE))

})

test_that("risib step works", {

    # autosome
    rf <- 0.01
    RF <- 4*rf/(1+6*rf)
    expect_equal(test_step("risib", 1, 1, rf, FALSE, FALSE, 0), log(1-RF))
    expect_equal(test_step("risib", 1, 2, rf, FALSE, FALSE, 0), log(RF))
    expect_equal(test_step("risib", 2, 1, rf, FALSE, FALSE, 0), log(RF))
    expect_equal(test_step("risib", 2, 2, rf, FALSE, FALSE, 0), log(1-RF))
    # X chr (see table 3 in Broman, Genetics 169:1133-1146, 2005)
    t11 <- (1+2*rf)/(1+4*rf)
    t12 <- 2*rf/(1+4*rf)
    t21 <- 4*rf/(1+4*rf)
    t22 <- 1/(1+4*rf)
    expect_equal(test_step("risib", 1, 1, rf, TRUE, TRUE, 0), log(t11))
    expect_equal(test_step("risib", 1, 2, rf, TRUE, TRUE, 0), log(t12))
    expect_equal(test_step("risib", 2, 1, rf, TRUE, TRUE, 0), log(t21))
    expect_equal(test_step("risib", 2, 2, rf, TRUE, TRUE, 0), log(t22))
    # X chr reverse
    expect_equal(test_step("risib", 1, 1, rf, TRUE, FALSE, 1), log(t22))
    expect_equal(test_step("risib", 1, 2, rf, TRUE, FALSE, 1), log(t21))
    expect_equal(test_step("risib", 2, 1, rf, TRUE, FALSE, 1), log(t12))
    expect_equal(test_step("risib", 2, 2, rf, TRUE, FALSE, 1), log(t11))

    # phase-known shouldn't matter
    rf <- 0.15
    RF <- 4*rf/(1+6*rf)
    expect_equal(test_step("risib", 1, 1, rf, FALSE, FALSE, 0, TRUE), log(1-RF))
    expect_equal(test_step("risib", 1, 2, rf, FALSE, FALSE, 0, TRUE), log(RF))
    expect_equal(test_step("risib", 2, 1, rf, FALSE, FALSE, 0, TRUE), log(RF))
    expect_equal(test_step("risib", 2, 2, rf, FALSE, FALSE, 0, TRUE), log(1-RF))
    # X chr
    t11 <- (1+2*rf)/(1+4*rf)
    t12 <- 2*rf/(1+4*rf)
    t21 <- 4*rf/(1+4*rf)
    t22 <- 1/(1+4*rf)
    expect_equal(test_step("risib", 1, 1, rf, TRUE, TRUE, 0, TRUE), log(t11))
    expect_equal(test_step("risib", 1, 2, rf, TRUE, TRUE, 0, TRUE), log(t12))
    expect_equal(test_step("risib", 2, 1, rf, TRUE, TRUE, 0, TRUE), log(t21))
    expect_equal(test_step("risib", 2, 2, rf, TRUE, TRUE, 0, TRUE), log(t22))
    # X chr reverse
    expect_equal(test_step("risib", 1, 1, rf, TRUE, FALSE, 1, TRUE), log(t22))
    expect_equal(test_step("risib", 1, 2, rf, TRUE, FALSE, 1, TRUE), log(t21))
    expect_equal(test_step("risib", 2, 1, rf, TRUE, FALSE, 1, TRUE), log(t12))
    expect_equal(test_step("risib", 2, 2, rf, TRUE, FALSE, 1, TRUE), log(t11))

    # errors
    expect_error(test_step("risib", 0, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("risib", 1, 0, rf, FALSE, FALSE, 0))
    expect_error(test_step("risib", 3, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("risib", 2, 3, rf, FALSE, FALSE, 0))
    # X chr
    expect_error(test_step("risib", 0, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("risib", 1, 0, rf, TRUE, TRUE, 0))
    expect_error(test_step("risib", 3, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("risib", 2, 3, rf, TRUE, TRUE, 0))
    # X chr reverse
    expect_error(test_step("risib", 0, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("risib", 1, 0, rf, TRUE, FALSE, 1))
    expect_error(test_step("risib", 3, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("risib", 2, 3, rf, TRUE, FALSE, 1))

    # phase-known shouldn't matter
    expect_error(test_step("risib", 0, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("risib", 1, 0, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("risib", 3, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("risib", 2, 3, rf, FALSE, FALSE, 0, TRUE))
    # X chr
    expect_error(test_step("risib", 0, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("risib", 1, 0, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("risib", 3, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("risib", 2, 3, rf, TRUE, TRUE, 0, TRUE))
    # X reverse
    expect_error(test_step("risib", 0, 1, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("risib", 1, 0, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("risib", 3, 1, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("risib", 2, 3, rf, TRUE, FALSE, 1, TRUE))

})

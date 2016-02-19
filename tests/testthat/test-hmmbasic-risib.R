context("basic HMM functions in risib")

test_that("risib check_geno works", {

    # autosome
    expect_true(test_check_geno("risib", 0, TRUE, FALSE, FALSE, 0))
    expect_true(test_check_geno("risib", 1, TRUE, FALSE, FALSE, 0))
    expect_true(test_check_geno("risib", 2, TRUE, FALSE, FALSE, 0))
    expect_true(test_check_geno("risib", 1, FALSE, FALSE, FALSE, 0))
    expect_true(test_check_geno("risib", 2, FALSE, FALSE, FALSE, 0))
    expect_false(test_check_geno("risib", 3, TRUE, FALSE, FALSE, 0))
    expect_false(test_check_geno("risib", 0, FALSE, FALSE, FALSE, 0))
    expect_false(test_check_geno("risib", 3, FALSE, FALSE, FALSE, 0))

    # X chromosome
    expect_true(test_check_geno("risib", 0, TRUE, TRUE, TRUE, 0))
    expect_true(test_check_geno("risib", 1, TRUE, TRUE, TRUE, 0))
    expect_true(test_check_geno("risib", 2, TRUE, TRUE, TRUE, 0))
    expect_true(test_check_geno("risib", 1, FALSE, TRUE, TRUE, 0))
    expect_true(test_check_geno("risib", 2, FALSE, TRUE, TRUE, 0))
    expect_false(test_check_geno("risib", 3, TRUE, TRUE, TRUE, 0))
    expect_false(test_check_geno("risib", 0, FALSE, TRUE, TRUE, 0))
    expect_false(test_check_geno("risib", 3, FALSE, TRUE, TRUE, 0))

    # X chr reverse
    expect_true(test_check_geno("risib", 0, TRUE, TRUE, TRUE, 1))
    expect_true(test_check_geno("risib", 1, TRUE, TRUE, TRUE, 1))
    expect_true(test_check_geno("risib", 2, TRUE, TRUE, TRUE, 1))
    expect_true(test_check_geno("risib", 1, FALSE, TRUE, TRUE, 1))
    expect_true(test_check_geno("risib", 2, FALSE, TRUE, TRUE, 1))
    expect_false(test_check_geno("risib", 3, TRUE, TRUE, TRUE, 1))
    expect_false(test_check_geno("risib", 0, FALSE, TRUE, TRUE, 1))
    expect_false(test_check_geno("risib", 3, FALSE, TRUE, TRUE, 1))

})

test_that("risib ngen works", {

    expect_equal(test_ngen("risib", FALSE), 2)
    expect_equal(test_ngen("risib", TRUE), 2)

})

test_that("risib possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("risib", FALSE, FALSE, 0), 1:2)
    # X chr
    expect_equal(test_possible_gen("risib", TRUE, TRUE, 0), 1:2)
    # X chr reverse
    expect_equal(test_possible_gen("risib", TRUE, TRUE, 1), 1:2)

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

})

test_that("risib emit works", {

    # autosome
    eps <- 0.01
    expect_equal(test_emit("risib", 0, 1, eps, integer(0), FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 0, 2, eps, integer(0), FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 1, 1, eps, integer(0), FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, integer(0), FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, integer(0), FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, integer(0), FALSE, FALSE, 0), log(1-eps))
    # X chr
    expect_equal(test_emit("risib", 0, 1, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 0, 2, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 1, 1, eps, integer(0), TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, integer(0), TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, integer(0), TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, integer(0), TRUE, TRUE, 0), log(1-eps))
    # X chr reverse
    expect_equal(test_emit("risib", 0, 1, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 0, 2, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 1, 1, eps, integer(0), TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("risib", 1, 2, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("risib", 2, 1, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("risib", 2, 2, eps, integer(0), TRUE, FALSE, 1), log(1-eps))

    # errors
    expect_equal(test_emit("risib", 0, 0, eps, integer(0), FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 0, 3, eps, integer(0), FALSE, FALSE, 0), 0)
    expect_equal(test_emit("risib", 3, 1, eps, integer(0), FALSE, FALSE, 0), 0)
    # X chr
    expect_equal(test_emit("risib", 0, 0, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 0, 3, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("risib", 3, 1, eps, integer(0), TRUE, TRUE, 0), 0)
    # X chr reverse
    expect_equal(test_emit("risib", 0, 0, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 0, 3, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("risib", 3, 1, eps, integer(0), TRUE, FALSE, 1), 0)

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

})

test_that("geno_names works", {
    expect_equal(geno_names("risib", c("B", "R"), FALSE), c("BB", "RR"))
    expect_equal(geno_names("risib", c("B", "R"), TRUE), c("BB", "RR"))
})

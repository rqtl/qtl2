
context("basic HMM functions in intercross")

test_that("intercross check_geno works", {

    # autosome
    expect_true(test_check_geno("f2", 0, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 2, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 3, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 4, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 4, TRUE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, FALSE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 2, FALSE, FALSE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 3, FALSE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 6, TRUE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 0, FALSE, FALSE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 4, FALSE, FALSE, FALSE, 0, FALSE))

    # X chromosome female, forward cross
    expect_true(test_check_geno("f2", 0, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 2, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 4, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 5, TRUE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, FALSE, TRUE, TRUE, 0, FALSE))
    expect_true(test_check_geno("f2", 2, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("f2", 3, TRUE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("f2", 0, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("f2", 3, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("f2", 4, FALSE, TRUE, TRUE, 0, FALSE))
    expect_error(test_check_geno("f2", 5, FALSE, TRUE, TRUE, 0, FALSE))

    # X chromosome female, reverse cross
    expect_true(test_check_geno("f2", 0, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 2, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 3, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 4, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 5, TRUE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 2, FALSE, TRUE, TRUE, 1, FALSE))
    expect_true(test_check_geno("f2", 3, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 1, TRUE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 6, TRUE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 0, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 1, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 4, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 5, FALSE, TRUE, TRUE, 1, FALSE))
    expect_error(test_check_geno("f2", 6, FALSE, TRUE, TRUE, 1, FALSE))

    # X chromosome male
    expect_true(test_check_geno("f2", 0, TRUE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, TRUE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 3, TRUE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 4, TRUE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 5, TRUE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 1, FALSE, TRUE, FALSE, 0, FALSE))
    expect_true(test_check_geno("f2", 3, FALSE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 2, TRUE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 6, TRUE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 0, FALSE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 2, FALSE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 4, FALSE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 5, FALSE, TRUE, FALSE, 0, FALSE))
    expect_error(test_check_geno("f2", 6, FALSE, TRUE, FALSE, 0, FALSE))

})

test_that("intercross all_geno works", {
    expect_equal(test_allgeno("f2", FALSE, FALSE), 1:3)
    expect_equal(test_allgeno("f2", TRUE, FALSE),  1:3)

    # phase-known
    expect_equal(test_allgeno("f2", FALSE, TRUE), 1:4)
    expect_equal(test_allgeno("f2", TRUE, TRUE),  1:4)
})

test_that("intercross geno works", {

    # autosome
    expect_equal(test_geno("f2", FALSE, FALSE, 0, FALSE), 1:3)
    # X female forward
    expect_equal(test_geno("f2", TRUE, TRUE, 0, FALSE), 1:2)
    # X female reverse
    expect_equal(test_geno("f2", TRUE, TRUE, 1, FALSE), 2:3)
    # X male
    expect_equal(test_geno("f2", TRUE, FALSE, 0, FALSE), c(1,3))
    # X male reverse
    expect_equal(test_geno("f2", TRUE, FALSE, 1, FALSE), c(1,3))

    # phase-known
    expect_equal(test_geno("f2", FALSE, FALSE, 0, TRUE), 1:4)
    expect_equal(test_geno("f2", TRUE, TRUE, 0, TRUE), 1:2)
    expect_equal(test_geno("f2", TRUE, TRUE, 1, TRUE), 3:4)
    expect_equal(test_geno("f2", TRUE, FALSE, 0, TRUE), c(1,4))
    expect_equal(test_geno("f2", TRUE, FALSE, 1, TRUE), c(1,4))
})


test_that("intercross nrec works", {

    # autosome
    expect_equal(test_nrec("f2", 1, 1, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2", 1, 2, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 1, 3, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 1, 4, FALSE, FALSE, 0), 2)
    expect_equal(test_nrec("f2", 2, 1, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 2, 2, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2", 2, 3, FALSE, FALSE, 0), 2)
    expect_equal(test_nrec("f2", 2, 4, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 3, 1, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 3, 2, FALSE, FALSE, 0), 2)
    expect_equal(test_nrec("f2", 3, 3, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2", 3, 4, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 4, 1, FALSE, FALSE, 0), 2)
    expect_equal(test_nrec("f2", 4, 2, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 4, 3, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 4, 4, FALSE, FALSE, 0), 0)
    # X female forward
    expect_equal(test_nrec("f2", 1, 1, TRUE, TRUE, 0), 0)
    expect_equal(test_nrec("f2", 1, 2, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("f2", 2, 1, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("f2", 2, 2, TRUE, TRUE, 0), 0)
    # X female reverse
    expect_equal(test_nrec("f2", 3, 3, TRUE, TRUE, 1), 0)
    expect_equal(test_nrec("f2", 3, 4, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("f2", 4, 3, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("f2", 4, 4, TRUE, TRUE, 1), 0)
    # X male
    expect_equal(test_nrec("f2", 1, 1, TRUE, FALSE, 0), 0)
    expect_equal(test_nrec("f2", 1, 4, TRUE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 4, 1, TRUE, FALSE, 0), 1)
    expect_equal(test_nrec("f2", 4, 4, TRUE, FALSE, 0), 0)
    # X male reverse
    expect_equal(test_nrec("f2", 1, 1, TRUE, FALSE, 1), 0)
    expect_equal(test_nrec("f2", 1, 4, TRUE, FALSE, 1), 1)
    expect_equal(test_nrec("f2", 4, 1, TRUE, FALSE, 1), 1)
    expect_equal(test_nrec("f2", 4, 4, TRUE, FALSE, 1), 0)

    # some errors
    # autosome
    expect_error(test_nrec("f2", 0, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("f2", 1, 5, FALSE, FALSE, 0))
    expect_error(test_nrec("f2", 5, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("f2", 1, 0, FALSE, FALSE, 0))
    # X female forward
    expect_error(test_nrec("f2", 0, 1, TRUE, TRUE, 0))
    expect_error(test_nrec("f2", 1, 3, TRUE, TRUE, 0))
    expect_error(test_nrec("f2", 3, 1, TRUE, TRUE, 0))
    expect_error(test_nrec("f2", 2, 0, TRUE, TRUE, 0))
    # X female reverse
    expect_error(test_nrec("f2", 0, 3, TRUE, TRUE, 1))
    expect_error(test_nrec("f2", 3, 1, TRUE, TRUE, 1))
    expect_error(test_nrec("f2", 2, 3, TRUE, TRUE, 1))
    expect_error(test_nrec("f2", 4, 0, TRUE, TRUE, 1))
    # X male
    expect_error(test_nrec("f2", 0, 1, TRUE, FALSE, 0))
    expect_error(test_nrec("f2", 2, 4, TRUE, FALSE, 0))
    expect_error(test_nrec("f2", 4, 3, TRUE, FALSE, 0))
    expect_error(test_nrec("f2", 4, 0, TRUE, FALSE, 0))
    # X male reverse
    expect_error(test_nrec("f2", 0, 1, TRUE, FALSE, 1))
    expect_error(test_nrec("f2", 1, 0, TRUE, FALSE, 1))
    expect_error(test_nrec("f2", 2, 1, TRUE, FALSE, 1))
    expect_error(test_nrec("f2", 4, 3, TRUE, FALSE, 1))

})

test_that("intercross init works", {

    # autosome
    expect_equal(test_init("f2", 1, FALSE, FALSE, 0), log(0.25))
    expect_equal(test_init("f2", 2, FALSE, FALSE, 0), log(0.5))
    expect_equal(test_init("f2", 3, FALSE, FALSE, 0), log(0.25))
    # X female forward
    for(i in 1:2)
        expect_equal(test_init("f2", i, TRUE, TRUE, 0), log(0.5))
    # X female reverse
    for(i in 2:3)
        expect_equal(test_init("f2", i, TRUE, TRUE, 1), log(0.5))
    # X male
    for(i in c(1,3))
        expect_equal(test_init("f2", i, TRUE, FALSE, 0), log(0.5))
    # X male reverse
    for(i in c(1,3))
        expect_equal(test_init("f2", i, TRUE, FALSE, 1), log(0.5))

    # phase-known
    # autosome
    for(i in 1:4)
        expect_equal(test_init("f2", i, FALSE, FALSE, 0, TRUE), log(0.25))
    # X female forward
    for(i in 1:2)
        expect_equal(test_init("f2", i, TRUE, TRUE, 0, TRUE), log(0.5))
    # X female reverse
    for(i in 3:4)
        expect_equal(test_init("f2", i, TRUE, TRUE, 1, TRUE), log(0.5))
    # X male
    for(i in c(1,4))
        expect_equal(test_init("f2", i, TRUE, FALSE, 0, TRUE), log(0.5))
    # X male reverse
    for(i in c(1,4))
        expect_equal(test_init("f2", i, TRUE, FALSE, 1, TRUE), log(0.5))

})

test_that("intercross emit works", {

    # autosome
    eps <- 0.01
    for(i in 1:3)
        expect_equal(test_emit("f2", 0, i, eps, FALSE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 1, 1, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 3, 1, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 4, 1, eps, FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2", 5, 1, eps, FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 2, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 2, 2, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 4, 2, eps, FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2", 5, 2, eps, FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2", 1, 3, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 2, 3, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 3, 3, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, FALSE, FALSE, 0), log(1-eps/2))
    # X female forward
    for(i in 1:2)
        expect_equal(test_emit("f2", 0, i, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, TRUE, 0), log(eps))
    expect_error(test_emit("f2", 3, 1, eps, TRUE, TRUE, 0))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 0), log(1-eps))
    expect_error(test_emit("f2", 3, 2, eps, TRUE, TRUE, 0))
    expect_equal(test_emit("f2", 4, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 2, eps, TRUE, TRUE, 0), 0)
    # X female reverse
    for(i in 2:3)
        expect_equal(test_emit("f2", 0, i, eps, TRUE, TRUE, 0), 0)
    expect_error(test_emit("f2", 1, 3, eps, TRUE, TRUE, 0))
    expect_equal(test_emit("f2", 2, 3, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 0), log(1-eps))
    expect_error(test_emit("f2", 3, 2, eps, TRUE, TRUE, 0))
    expect_equal(test_emit("f2", 4, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 2, eps, TRUE, TRUE, 0), 0)
    # X male
    for(i in c(1,3))
        expect_equal(test_emit("f2", 0, i, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_error(test_emit("f2", 2, 1, eps, TRUE, FALSE, 0))
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, FALSE, 0), log(eps))
    expect_error(test_emit("f2", 2, 3, eps, TRUE, FALSE, 0))
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, FALSE, 0), log(1-eps))

######################################################################
# finished here, but still some problems above
#
# for X chromosome in intercross, I think we could allow any observed
# genotype for either sex or direction, but treat invalid ones as NA
# (no need to throw error)
######################################################################

    # phase-known
    eps <- 0.001
    expect_equal(test_emit("f2", 0, 1, eps, FALSE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 0, 2, eps, FALSE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 1, 1, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 1, 2, eps, FALSE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 1, eps, FALSE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    # X female
    expect_equal(test_emit("f2", 0, 1, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 0, 2, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    # X male
    expect_equal(test_emit("f2", 0, 1, eps, TRUE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 0, 3, eps, TRUE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, FALSE, 0, TRUE), log(1-eps))

    # errors
    expect_error(test_emit("f2", 0, 0, eps, FALSE, FALSE, 0))
    expect_error(test_emit("f2", 0, 3, eps, FALSE, FALSE, 0))
    expect_error(test_emit("f2", 3, 1, eps, FALSE, FALSE, 0))
    # X female
    expect_error(test_emit("f2", 0, 0, eps, TRUE, TRUE, 0))
    expect_error(test_emit("f2", 0, 3, eps, TRUE, TRUE, 0))
    expect_error(test_emit("f2", 3, 1, eps, TRUE, TRUE, 0))
    # X male
    expect_error(test_emit("f2", 0, 0, eps, TRUE, FALSE, 0))
    expect_error(test_emit("f2", 0, 2, eps, TRUE, FALSE, 0))
    expect_error(test_emit("f2", 2, 1, eps, TRUE, FALSE, 0))

    # errors
    expect_error(test_emit("f2", 0, 0, eps, FALSE, FALSE, 0, TRUE))
    expect_error(test_emit("f2", 0, 3, eps, FALSE, FALSE, 0, TRUE))
    expect_error(test_emit("f2", 3, 1, eps, FALSE, FALSE, 0, TRUE))
    # X female
    expect_error(test_emit("f2", 0, 0, eps, TRUE, TRUE, 0, TRUE))
    expect_error(test_emit("f2", 0, 3, eps, TRUE, TRUE, 0, TRUE))
    expect_error(test_emit("f2", 3, 1, eps, TRUE, TRUE, 0, TRUE))
    # X male
    expect_error(test_emit("f2", 0, 0, eps, TRUE, FALSE, 0, TRUE))
    expect_error(test_emit("f2", 0, 2, eps, TRUE, FALSE, 0, TRUE))
    expect_error(test_emit("f2", 2, 1, eps, TRUE, FALSE, 0, TRUE))

})

test_that("intercross step works", {

    # autosome
    rf <- 0.01
    expect_equal(test_step("f2", 1, 1, rf, FALSE, FALSE, 0), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, FALSE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, FALSE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, FALSE, FALSE, 0), log(1-rf))
    # X female
    expect_equal(test_step("f2", 1, 1, rf, TRUE, TRUE, 0), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, TRUE, TRUE, 0), log(1-rf))
    # X male
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 0), log(1-rf))
    expect_equal(test_step("f2", 1, 3, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 3, 1, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 3, 3, rf, TRUE, FALSE, 0), log(1-rf))

    # phase-known
    rf <- 0.15
    expect_equal(test_step("f2", 1, 1, rf, FALSE, FALSE, 0, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, FALSE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, FALSE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, FALSE, FALSE, 0, TRUE), log(1-rf))
    # X female
    expect_equal(test_step("f2", 1, 1, rf, TRUE, TRUE, 0, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, TRUE, TRUE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, TRUE, TRUE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, TRUE, TRUE, 0, TRUE), log(1-rf))
    # X male
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 0, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 3, rf, TRUE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 3, 1, rf, TRUE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 3, 3, rf, TRUE, FALSE, 0, TRUE), log(1-rf))

    # errors
    expect_error(test_step("f2", 0, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 1, 0, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 3, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 2, 3, rf, FALSE, FALSE, 0))
    # X female
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 3, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 2, 3, rf, TRUE, TRUE, 0))
    # X male
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 3, 2, rf, TRUE, FALSE, 0))

    # phase-known
    expect_error(test_step("f2", 0, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 3, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 2, 3, rf, FALSE, FALSE, 0, TRUE))
    # X female
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 3, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 2, 3, rf, TRUE, TRUE, 0, TRUE))
    # X male
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 3, 2, rf, TRUE, FALSE, 0, TRUE))

})

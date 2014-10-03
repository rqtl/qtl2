
context("basic HMM functions in intercross")

test_that("intercross check_geno works", {

    # autosome
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0, FALSE))
    for(i in 1:3)
        expect_true(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0, FALSE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0, FALSE))
    for(i in c(0, 4))
        expect_error(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0, FALSE))

    # X chromosome female, forward cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0, FALSE))
    for(i in 1:2)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0, FALSE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0, FALSE))
    for(i in c(0, 3, 4))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0, FALSE))

    # X chromosome female, reverse cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1, FALSE))
    for(i in 2:3)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1, FALSE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1, FALSE))
    for(i in c(0, 1, 4))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1, FALSE))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0, FALSE))
    for(i in c(1,3))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0, FALSE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0, FALSE))
    for(i in c(0, 2, 4))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0, FALSE))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1, FALSE))
    for(i in c(1,3))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1, FALSE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1, FALSE))
    for(i in c(0, 2, 4))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1, FALSE))

})

test_that("intercross check_geno works in phase-known case", {

    # autosome
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0, TRUE))
    for(i in 1:4)
        expect_true(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0, TRUE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0, TRUE))
    for(i in c(0, 5))
        expect_error(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0, TRUE))

    # X chromosome female, forward cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0, TRUE))
    for(i in 1:2)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0, TRUE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0, TRUE))
    for(i in c(0, 3, 4))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0, TRUE))

    # X chromosome female, reverse cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1, TRUE))
    for(i in 3:4)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1, TRUE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1, TRUE))
    for(i in c(0, 1, 2, 5))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1, TRUE))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0, TRUE))
    for(i in c(1,4))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0, TRUE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0, TRUE))
    for(i in c(0, 2, 3, 5))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0, TRUE))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1, TRUE))
    for(i in c(1,4))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1, TRUE))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1, TRUE))
    for(i in c(0, 2, 3, 5))
        expect_error(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1, TRUE))

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
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, TRUE, 0), log(eps))

    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 4, 2, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 2, eps, TRUE, TRUE, 0), log(1-eps))

    # X female reverse
    for(i in 2:3)
        expect_equal(test_emit("f2", 0, i, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 2, 3, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 4, 2, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 5, 2, eps, TRUE, TRUE, 1), 0)

    # X male
    for(i in c(1,3))
        expect_equal(test_emit("f2", 0, i, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 3, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, FALSE, 0), log(1-eps))

    # X male reverse
    for(i in c(1,3))
        expect_equal(test_emit("f2", 0, i, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 2, 3, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, FALSE, 1), log(1-eps))

    # phase-known
    # autosome
    eps <- 0.02
    for(i in 1:4)
        expect_equal(test_emit("f2", 0, i, eps, FALSE, FALSE, 0, TRUE), 0)

    expect_equal(test_emit("f2", 1, 1, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
    expect_equal(test_emit("f2", 3, 1, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
    expect_equal(test_emit("f2", 4, 1, eps, FALSE, FALSE, 0, TRUE), log(1-eps/2))
    expect_equal(test_emit("f2", 5, 1, eps, FALSE, FALSE, 0, TRUE), log(eps))

    for(i in 2:3) {
        expect_equal(test_emit("f2", 1, i, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
        expect_equal(test_emit("f2", 2, i, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
        expect_equal(test_emit("f2", 3, i, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
        expect_equal(test_emit("f2", 4, i, eps, FALSE, FALSE, 0, TRUE), log(1-eps/2))
        expect_equal(test_emit("f2", 5, i, eps, FALSE, FALSE, 0, TRUE), log(1-eps/2))
    }

    expect_equal(test_emit("f2", 1, 4, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
    expect_equal(test_emit("f2", 2, 4, eps, FALSE, FALSE, 0, TRUE), log(eps/2))
    expect_equal(test_emit("f2", 3, 4, eps, FALSE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 4, 4, eps, FALSE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 5, 4, eps, FALSE, FALSE, 0, TRUE), log(1-eps/2))

    # X female forward
    for(i in 1:2)
        expect_equal(test_emit("f2", 0, i, eps, TRUE, TRUE, 0, TRUE), 0)

    expect_equal(test_emit("f2", 1, 1, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, TRUE, 0, TRUE), log(eps))

    expect_equal(test_emit("f2", 1, 2, eps, TRUE, TRUE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, TRUE, TRUE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 4, 2, eps, TRUE, TRUE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 5, 2, eps, TRUE, TRUE, 0, TRUE), log(1-eps))

    # X female reverse
    for(i in 3:4)
        expect_equal(test_emit("f2", 0, i, eps, TRUE, TRUE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 1, 4, eps, TRUE, TRUE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 2, 4, eps, TRUE, TRUE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 3, 4, eps, TRUE, TRUE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 4, 4, eps, TRUE, TRUE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 5, 4, eps, TRUE, TRUE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 1, 3, eps, TRUE, TRUE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 2, 3, eps, TRUE, TRUE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 3, 3, eps, TRUE, TRUE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 4, 3, eps, TRUE, TRUE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 5, 3, eps, TRUE, TRUE, 1, TRUE), 0)

    # X male
    for(i in c(1,4))
        expect_equal(test_emit("f2", 0, i, eps, TRUE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 1, 4, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 4, eps, TRUE, FALSE, 0, TRUE), 0)
    expect_equal(test_emit("f2", 3, 4, eps, TRUE, FALSE, 0, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 4, 4, eps, TRUE, FALSE, 0, TRUE), log(eps))
    expect_equal(test_emit("f2", 5, 4, eps, TRUE, FALSE, 0, TRUE), log(1-eps))

    # X male reverse
    for(i in c(1,4))
        expect_equal(test_emit("f2", 0, i, eps, TRUE, FALSE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 1, 1, eps, TRUE, FALSE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, TRUE, FALSE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 3, 1, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 4, 1, eps, TRUE, FALSE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 5, 1, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 1, 4, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 2, 4, eps, TRUE, FALSE, 1, TRUE), 0)
    expect_equal(test_emit("f2", 3, 4, eps, TRUE, FALSE, 1, TRUE), log(1-eps))
    expect_equal(test_emit("f2", 4, 4, eps, TRUE, FALSE, 1, TRUE), log(eps))
    expect_equal(test_emit("f2", 5, 4, eps, TRUE, FALSE, 1, TRUE), log(1-eps))

})

test_that("intercross step works", {

    # autosome
    rf <- 0.01
    expect_equal(test_step("f2", 1, 1, rf, FALSE, FALSE, 0), log((1-rf)^2))
    expect_equal(test_step("f2", 1, 2, rf, FALSE, FALSE, 0), log(2*rf*(1-rf)))
    expect_equal(test_step("f2", 1, 3, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2", 2, 1, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2", 2, 2, rf, FALSE, FALSE, 0), log((1-rf)^2+rf^2))
    expect_equal(test_step("f2", 2, 3, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2", 3, 1, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2", 3, 2, rf, FALSE, FALSE, 0), log(2*rf*(1-rf)))
    expect_equal(test_step("f2", 3, 3, rf, FALSE, FALSE, 0), log((1-rf)^2))
    # X female
    expect_equal(test_step("f2", 1, 1, rf, TRUE, TRUE, 0), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, TRUE, TRUE, 0), log(1-rf))
    # X female reverse
    expect_equal(test_step("f2", 3, 3, rf, TRUE, TRUE, 1), log(1-rf))
    expect_equal(test_step("f2", 3, 2, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2", 2, 3, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, TRUE, TRUE, 1), log(1-rf))
    # X male
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 0), log(1-rf))
    expect_equal(test_step("f2", 1, 3, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 3, 1, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 3, 3, rf, TRUE, FALSE, 0), log(1-rf))
    # X male reverse
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 1), log(1-rf))
    expect_equal(test_step("f2", 1, 3, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2", 3, 1, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2", 3, 3, rf, TRUE, FALSE, 1), log(1-rf))

    # phase-known
    rf <- 0.15
    expect_equal(test_step("f2", 1, 1, rf, FALSE, FALSE, 0, TRUE), log((1-rf)^2))
    expect_equal(test_step("f2", 1, 2, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 1, 3, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 1, 4, rf, FALSE, FALSE, 0, TRUE), log(rf^2))
    expect_equal(test_step("f2", 2, 1, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 2, 2, rf, FALSE, FALSE, 0, TRUE), log((1-rf)^2))
    expect_equal(test_step("f2", 2, 3, rf, FALSE, FALSE, 0, TRUE), log(rf^2))
    expect_equal(test_step("f2", 2, 4, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 3, 1, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 3, 2, rf, FALSE, FALSE, 0, TRUE), log(rf^2))
    expect_equal(test_step("f2", 3, 3, rf, FALSE, FALSE, 0, TRUE), log((1-rf)^2))
    expect_equal(test_step("f2", 3, 4, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 4, 1, rf, FALSE, FALSE, 0, TRUE), log(rf^2))
    expect_equal(test_step("f2", 4, 2, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 4, 3, rf, FALSE, FALSE, 0, TRUE), log(rf*(1-rf)))
    expect_equal(test_step("f2", 4, 4, rf, FALSE, FALSE, 0, TRUE), log((1-rf)^2))
    # X female
    expect_equal(test_step("f2", 1, 1, rf, TRUE, TRUE, 0, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 2, rf, TRUE, TRUE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 1, rf, TRUE, TRUE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 2, 2, rf, TRUE, TRUE, 0, TRUE), log(1-rf))
    # X female reverse
    expect_equal(test_step("f2", 3, 3, rf, TRUE, TRUE, 1, TRUE), log(1-rf))
    expect_equal(test_step("f2", 3, 4, rf, TRUE, TRUE, 1, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 3, rf, TRUE, TRUE, 1, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 4, rf, TRUE, TRUE, 1, TRUE), log(1-rf))
    # X male
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 0, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 4, rf, TRUE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 1, rf, TRUE, FALSE, 0, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 4, rf, TRUE, FALSE, 0, TRUE), log(1-rf))
    # X male reverse
    expect_equal(test_step("f2", 1, 1, rf, TRUE, FALSE, 1, TRUE), log(1-rf))
    expect_equal(test_step("f2", 1, 4, rf, TRUE, FALSE, 1, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 1, rf, TRUE, FALSE, 1, TRUE), log(rf))
    expect_equal(test_step("f2", 4, 4, rf, TRUE, FALSE, 1, TRUE), log(1-rf))

    # errors
    expect_error(test_step("f2", 0, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 1, 0, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 4, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2", 2, 4, rf, FALSE, FALSE, 0))
    # X female
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 3, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2", 2, 3, rf, TRUE, TRUE, 0))
    # X female reverse
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2", 1, 1, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2", 2, 1, rf, TRUE, TRUE, 1))
    # X male
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2", 4, 2, rf, TRUE, FALSE, 0))
    # X male reverse
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2", 4, 2, rf, TRUE, FALSE, 1))

    # errors
    expect_error(test_step("f2", 0, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 5, 1, rf, FALSE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 2, 5, rf, FALSE, FALSE, 0, TRUE))
    # X female
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 3, 1, rf, TRUE, TRUE, 0, TRUE))
    expect_error(test_step("f2", 2, 3, rf, TRUE, TRUE, 0, TRUE))
    # X female reverse
    expect_error(test_step("f2", 0, 1, rf, TRUE, TRUE, 1, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, TRUE, 1, TRUE))
    expect_error(test_step("f2", 1, 1, rf, TRUE, TRUE, 1, TRUE))
    expect_error(test_step("f2", 2, 2, rf, TRUE, TRUE, 1, TRUE))
    # X male
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 0, TRUE))
    expect_error(test_step("f2", 4, 2, rf, TRUE, FALSE, 0, TRUE))
    # X male reverse
    expect_error(test_step("f2", 0, 1, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("f2", 1, 0, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("f2", 2, 1, rf, TRUE, FALSE, 1, TRUE))
    expect_error(test_step("f2", 4, 2, rf, TRUE, FALSE, 1, TRUE))

})

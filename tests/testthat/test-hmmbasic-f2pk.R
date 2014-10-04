
context("basic HMM functions in phase-known intercross")

test_that("p-k intercross check_geno works", {

    # autosome
    for(i in 0:5)
        expect_true(test_check_geno("f2pk", i, TRUE, FALSE, FALSE, 0))
    for(i in 1:4)
        expect_true(test_check_geno("f2pk", i, FALSE, FALSE, FALSE, 0))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2pk", i, TRUE, FALSE, FALSE, 0))
    for(i in c(0, 5))
        expect_error(test_check_geno("f2pk", i, FALSE, FALSE, FALSE, 0))

    # X chromosome female, forward cross
    for(i in 0:5)
        expect_true(test_check_geno("f2pk", i, TRUE, TRUE, TRUE, 0))
    for(i in 1:2)
        expect_true(test_check_geno("f2pk", i, FALSE, TRUE, TRUE, 0))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2pk", i, TRUE, TRUE, TRUE, 0))
    for(i in c(0, 3, 4))
        expect_error(test_check_geno("f2pk", i, FALSE, TRUE, TRUE, 0))

    # X chromosome female, reverse cross
    for(i in 0:5)
        expect_true(test_check_geno("f2pk", i, TRUE, TRUE, TRUE, 1))
    for(i in 3:4)
        expect_true(test_check_geno("f2pk", i, FALSE, TRUE, TRUE, 1))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2pk", i, TRUE, TRUE, TRUE, 1))
    for(i in c(0, 1, 2, 5))
        expect_error(test_check_geno("f2pk", i, FALSE, TRUE, TRUE, 1))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2pk", i, TRUE, TRUE, FALSE, 0))
    for(i in c(1,4))
        expect_true(test_check_geno("f2pk", i, FALSE, TRUE, FALSE, 0))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2pk", i, TRUE, TRUE, FALSE, 0))
    for(i in c(0, 2, 3, 5))
        expect_error(test_check_geno("f2pk", i, FALSE, TRUE, FALSE, 0))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2pk", i, TRUE, TRUE, FALSE, 1))
    for(i in c(1,4))
        expect_true(test_check_geno("f2pk", i, FALSE, TRUE, FALSE, 1))
    for(i in c(-1, 6))
        expect_error(test_check_geno("f2pk", i, TRUE, TRUE, FALSE, 1))
    for(i in c(0, 2, 3, 5))
        expect_error(test_check_geno("f2pk", i, FALSE, TRUE, FALSE, 1))

})


test_that("p-k intercross n_geno works", {

    expect_equal(test_ngen("f2pk", FALSE), 4)
    expect_equal(test_ngen("f2pk", TRUE),  4)
})

test_that("p-k intercross geno_index works", {

    expect_equal(test_possible_gen("f2pk", FALSE, FALSE, 0), 1:4)
    expect_equal(test_possible_gen("f2pk", TRUE, TRUE, 0), 1:2)
    expect_equal(test_possible_gen("f2pk", TRUE, TRUE, 1), 3:4)
    expect_equal(test_possible_gen("f2pk", TRUE, FALSE, 0), c(1,4))
    expect_equal(test_possible_gen("f2pk", TRUE, FALSE, 1), c(1,4))
})


test_that("p-k intercross nrec works", {

    # autosome
    expect_equal(test_nrec("f2pk", 1, 1, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2pk", 1, 2, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 1, 3, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 1, 4, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 2, 1, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 2, 2, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2pk", 2, 3, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 2, 4, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 3, 1, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 3, 2, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 3, 3, FALSE, FALSE, 0), 0)
    expect_equal(test_nrec("f2pk", 3, 4, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 4, 1, FALSE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 4, 2, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 4, 3, FALSE, FALSE, 0), 1/2)
    expect_equal(test_nrec("f2pk", 4, 4, FALSE, FALSE, 0), 0)
    # X female forward
    expect_equal(test_nrec("f2pk", 1, 1, TRUE, TRUE, 0), 0)
    expect_equal(test_nrec("f2pk", 1, 2, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("f2pk", 2, 1, TRUE, TRUE, 0), 1)
    expect_equal(test_nrec("f2pk", 2, 2, TRUE, TRUE, 0), 0)
    # X female reverse
    expect_equal(test_nrec("f2pk", 3, 3, TRUE, TRUE, 1), 0)
    expect_equal(test_nrec("f2pk", 3, 4, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("f2pk", 4, 3, TRUE, TRUE, 1), 1)
    expect_equal(test_nrec("f2pk", 4, 4, TRUE, TRUE, 1), 0)
    # X male
    expect_equal(test_nrec("f2pk", 1, 1, TRUE, FALSE, 0), 0)
    expect_equal(test_nrec("f2pk", 1, 4, TRUE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 4, 1, TRUE, FALSE, 0), 1)
    expect_equal(test_nrec("f2pk", 4, 4, TRUE, FALSE, 0), 0)
    # X male reverse
    expect_equal(test_nrec("f2pk", 1, 1, TRUE, FALSE, 1), 0)
    expect_equal(test_nrec("f2pk", 1, 4, TRUE, FALSE, 1), 1)
    expect_equal(test_nrec("f2pk", 4, 1, TRUE, FALSE, 1), 1)
    expect_equal(test_nrec("f2pk", 4, 4, TRUE, FALSE, 1), 0)

    # some errors
    # autosome
    expect_error(test_nrec("f2pk", 0, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("f2pk", 1, 5, FALSE, FALSE, 0))
    expect_error(test_nrec("f2pk", 5, 1, FALSE, FALSE, 0))
    expect_error(test_nrec("f2pk", 1, 0, FALSE, FALSE, 0))
    # X female forward
    expect_error(test_nrec("f2pk", 0, 1, TRUE, TRUE, 0))
    expect_error(test_nrec("f2pk", 1, 3, TRUE, TRUE, 0))
    expect_error(test_nrec("f2pk", 3, 1, TRUE, TRUE, 0))
    expect_error(test_nrec("f2pk", 2, 0, TRUE, TRUE, 0))
    # X female reverse
    expect_error(test_nrec("f2pk", 0, 3, TRUE, TRUE, 1))
    expect_error(test_nrec("f2pk", 3, 1, TRUE, TRUE, 1))
    expect_error(test_nrec("f2pk", 2, 3, TRUE, TRUE, 1))
    expect_error(test_nrec("f2pk", 4, 0, TRUE, TRUE, 1))
    # X male
    expect_error(test_nrec("f2pk", 0, 1, TRUE, FALSE, 0))
    expect_error(test_nrec("f2pk", 2, 4, TRUE, FALSE, 0))
    expect_error(test_nrec("f2pk", 4, 3, TRUE, FALSE, 0))
    expect_error(test_nrec("f2pk", 4, 0, TRUE, FALSE, 0))
    # X male reverse
    expect_error(test_nrec("f2pk", 0, 1, TRUE, FALSE, 1))
    expect_error(test_nrec("f2pk", 1, 0, TRUE, FALSE, 1))
    expect_error(test_nrec("f2pk", 2, 1, TRUE, FALSE, 1))
    expect_error(test_nrec("f2pk", 4, 3, TRUE, FALSE, 1))

})

test_that("p-k intercross init works", {

    # autosome
    for(i in 1:4)
        expect_equal(test_init("f2pk", i, FALSE, FALSE, 0), log(0.25))
    # X female forward
    for(i in 1:2)
        expect_equal(test_init("f2pk", i, TRUE, TRUE, 0), log(0.5))
    # X female reverse
    for(i in 3:4)
        expect_equal(test_init("f2pk", i, TRUE, TRUE, 1), log(0.5))
    # X male
    for(i in c(1,4))
        expect_equal(test_init("f2pk", i, TRUE, FALSE, 0), log(0.5))
    # X male reverse
    for(i in c(1,4))
        expect_equal(test_init("f2pk", i, TRUE, FALSE, 1), log(0.5))

})

test_that("p-k intercross emit works", {

    # autosome
    eps <- 0.02
    for(i in 1:4)
        expect_equal(test_emit("f2pk", 0, i, eps, FALSE, FALSE, 0), 0)

    expect_equal(test_emit("f2pk", 1, 1, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 2, 1, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2pk", 3, 1, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2pk", 4, 1, eps, FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2pk", 5, 1, eps, FALSE, FALSE, 0), log(eps))

    for(i in 2:3) {
        expect_equal(test_emit("f2pk", 1, i, eps, FALSE, FALSE, 0), log(eps/2))
        expect_equal(test_emit("f2pk", 2, i, eps, FALSE, FALSE, 0), log(1-eps))
        expect_equal(test_emit("f2pk", 3, i, eps, FALSE, FALSE, 0), log(eps/2))
        expect_equal(test_emit("f2pk", 4, i, eps, FALSE, FALSE, 0), log(1-eps/2))
        expect_equal(test_emit("f2pk", 5, i, eps, FALSE, FALSE, 0), log(1-eps/2))
    }

    expect_equal(test_emit("f2pk", 1, 4, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2pk", 2, 4, eps, FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2pk", 3, 4, eps, FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 4, 4, eps, FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2pk", 5, 4, eps, FALSE, FALSE, 0), log(1-eps/2))

    # X female forward
    for(i in 1:2)
        expect_equal(test_emit("f2pk", 0, i, eps, TRUE, TRUE, 0), 0)

    expect_equal(test_emit("f2pk", 1, 1, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 2, 1, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2pk", 3, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2pk", 4, 1, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2pk", 5, 1, eps, TRUE, TRUE, 0), log(eps))

    expect_equal(test_emit("f2pk", 1, 2, eps, TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2pk", 2, 2, eps, TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 3, 2, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2pk", 4, 2, eps, TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2pk", 5, 2, eps, TRUE, TRUE, 0), log(1-eps))

    # X female reverse
    for(i in 3:4)
        expect_equal(test_emit("f2pk", 0, i, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2pk", 1, 4, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2pk", 2, 4, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2pk", 3, 4, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 4, 4, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2pk", 5, 4, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2pk", 1, 3, eps, TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2pk", 2, 3, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 3, 3, eps, TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2pk", 4, 3, eps, TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 5, 3, eps, TRUE, TRUE, 1), 0)

    # X male
    for(i in c(1,4))
        expect_equal(test_emit("f2pk", 0, i, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2pk", 1, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 2, 1, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2pk", 3, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2pk", 4, 1, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 5, 1, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2pk", 1, 4, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2pk", 2, 4, eps, TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2pk", 3, 4, eps, TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2pk", 4, 4, eps, TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2pk", 5, 4, eps, TRUE, FALSE, 0), log(1-eps))

    # X male reverse
    for(i in c(1,4))
        expect_equal(test_emit("f2pk", 0, i, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2pk", 1, 1, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 2, 1, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2pk", 3, 1, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2pk", 4, 1, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 5, 1, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2pk", 1, 4, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2pk", 2, 4, eps, TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2pk", 3, 4, eps, TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2pk", 4, 4, eps, TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2pk", 5, 4, eps, TRUE, FALSE, 1), log(1-eps))

})

test_that("p-k intercross step works", {

    # autosome
    rf <- 0.15
    expect_equal(test_step("f2pk", 1, 1, rf, FALSE, FALSE, 0), log((1-rf)^2))
    expect_equal(test_step("f2pk", 1, 2, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 1, 3, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 1, 4, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2pk", 2, 1, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 2, 2, rf, FALSE, FALSE, 0), log((1-rf)^2))
    expect_equal(test_step("f2pk", 2, 3, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2pk", 2, 4, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 3, 1, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 3, 2, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2pk", 3, 3, rf, FALSE, FALSE, 0), log((1-rf)^2))
    expect_equal(test_step("f2pk", 3, 4, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 4, 1, rf, FALSE, FALSE, 0), log(rf^2))
    expect_equal(test_step("f2pk", 4, 2, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 4, 3, rf, FALSE, FALSE, 0), log(rf*(1-rf)))
    expect_equal(test_step("f2pk", 4, 4, rf, FALSE, FALSE, 0), log((1-rf)^2))
    # X female
    expect_equal(test_step("f2pk", 1, 1, rf, TRUE, TRUE, 0), log(1-rf))
    expect_equal(test_step("f2pk", 1, 2, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2pk", 2, 1, rf, TRUE, TRUE, 0), log(rf))
    expect_equal(test_step("f2pk", 2, 2, rf, TRUE, TRUE, 0), log(1-rf))
    # X female reverse
    expect_equal(test_step("f2pk", 3, 3, rf, TRUE, TRUE, 1), log(1-rf))
    expect_equal(test_step("f2pk", 3, 4, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2pk", 4, 3, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2pk", 4, 4, rf, TRUE, TRUE, 1), log(1-rf))
    # X male
    expect_equal(test_step("f2pk", 1, 1, rf, TRUE, FALSE, 0), log(1-rf))
    expect_equal(test_step("f2pk", 1, 4, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2pk", 4, 1, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2pk", 4, 4, rf, TRUE, FALSE, 0), log(1-rf))
    # X male reverse
    expect_equal(test_step("f2pk", 1, 1, rf, TRUE, FALSE, 1), log(1-rf))
    expect_equal(test_step("f2pk", 1, 4, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2pk", 4, 1, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2pk", 4, 4, rf, TRUE, FALSE, 1), log(1-rf))

    # errors
    expect_error(test_step("f2pk", 0, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2pk", 1, 0, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2pk", 5, 1, rf, FALSE, FALSE, 0))
    expect_error(test_step("f2pk", 2, 5, rf, FALSE, FALSE, 0))
    # X female
    expect_error(test_step("f2pk", 0, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2pk", 1, 0, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2pk", 3, 1, rf, TRUE, TRUE, 0))
    expect_error(test_step("f2pk", 2, 3, rf, TRUE, TRUE, 0))
    # X female reverse
    expect_error(test_step("f2pk", 0, 1, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2pk", 1, 0, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2pk", 1, 1, rf, TRUE, TRUE, 1))
    expect_error(test_step("f2pk", 2, 2, rf, TRUE, TRUE, 1))
    # X male
    expect_error(test_step("f2pk", 0, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2pk", 1, 0, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2pk", 2, 1, rf, TRUE, FALSE, 0))
    expect_error(test_step("f2pk", 4, 2, rf, TRUE, FALSE, 0))
    # X male reverse
    expect_error(test_step("f2pk", 0, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2pk", 1, 0, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2pk", 2, 1, rf, TRUE, FALSE, 1))
    expect_error(test_step("f2pk", 4, 2, rf, TRUE, FALSE, 1))

})

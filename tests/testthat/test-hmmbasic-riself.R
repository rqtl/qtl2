
context("basic HMM functions in riself")

test_that("riself check_geno works", {

    expect_true(test_check_geno("riself", 0, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("riself", 1, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("riself", 2, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("riself", 1, FALSE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("riself", 2, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("riself", 3, TRUE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("riself", 0, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("riself", 3, FALSE, FALSE, FALSE, integer(0)))

})

test_that("riself ngen works", {
    expect_equal(test_ngen("riself", FALSE), 2)
})

test_that("riself possible_gen works", {

    expect_equal(test_possible_gen("riself", FALSE, FALSE, integer(0)), 1:2)

})

test_that("riself nrec works", {

    expect_equal(test_nrec("riself", 1, 1, FALSE, FALSE, integer(0)), 0)
    expect_equal(test_nrec("riself", 1, 2, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("riself", 2, 1, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("riself", 2, 2, FALSE, FALSE, integer(0)), 0)

})

test_that("riself init works", {

    expect_equal(test_init("riself", 1, FALSE, FALSE, integer(0)), log(0.5))
    expect_equal(test_init("riself", 2, FALSE, FALSE, integer(0)), log(0.5))

})

test_that("riself emit works", {

    eps <- 0.01
    expect_equal(test_emit("riself", 0, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("riself", 0, 2, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("riself", 1, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))
    expect_equal(test_emit("riself", 1, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("riself", 2, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("riself", 2, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))

    # errors
    expect_equal(test_emit("riself", 0, 0, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("riself", 0, 3, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("riself", 3, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)

})

test_that("riself step works", {

    rf <- 0.01
    RF <- 2*rf/(1+2*rf)
    expect_equal(test_step("riself", 1, 1, rf, FALSE, FALSE, integer(0)), log(1-RF))
    expect_equal(test_step("riself", 1, 2, rf, FALSE, FALSE, integer(0)), log(RF))
    expect_equal(test_step("riself", 2, 1, rf, FALSE, FALSE, integer(0)), log(RF))
    expect_equal(test_step("riself", 2, 2, rf, FALSE, FALSE, integer(0)), log(1-RF))

})

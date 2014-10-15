
context("basic HMM functions in backcross")

test_that("backcross check_geno works", {

    # autosome
    expect_true(test_check_geno("bc", 0, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 1, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 2, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 1, FALSE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 2, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 3, TRUE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 0, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 3, FALSE, FALSE, FALSE, integer(0)))

    # X chromosome female
    expect_true(test_check_geno("bc", 0, TRUE, TRUE, TRUE, integer(0)))
    expect_true(test_check_geno("bc", 1, TRUE, TRUE, TRUE, integer(0)))
    expect_true(test_check_geno("bc", 2, TRUE, TRUE, TRUE, integer(0)))
    expect_true(test_check_geno("bc", 1, FALSE, TRUE, TRUE, integer(0)))
    expect_true(test_check_geno("bc", 2, FALSE, TRUE, TRUE, integer(0)))
    expect_false(test_check_geno("bc", 3, TRUE, TRUE, TRUE, integer(0)))
    expect_false(test_check_geno("bc", 0, FALSE, TRUE, TRUE, integer(0)))
    expect_false(test_check_geno("bc", 3, FALSE, TRUE, TRUE, integer(0)))

    # X chromosome male
    expect_true(test_check_geno("bc", 0, TRUE, TRUE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 1, TRUE, TRUE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 3, TRUE, TRUE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 3, FALSE, TRUE, FALSE, integer(0)))
    expect_true(test_check_geno("bc", 4, FALSE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 2, TRUE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 4, TRUE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 0, FALSE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 1, FALSE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 2, FALSE, TRUE, FALSE, integer(0)))
    expect_false(test_check_geno("bc", 5, FALSE, TRUE, FALSE, integer(0)))

})

test_that("backcross n_gen works", {

    expect_equal(test_ngen("bc", FALSE), 2)
    expect_equal(test_ngen("bc", TRUE), 4)

})

test_that("backcross possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("bc", FALSE, FALSE, integer(0)), 1:2)
    # X female
    expect_equal(test_possible_gen("bc", TRUE, TRUE, integer(0)), 1:2)
    # X male
    expect_equal(test_possible_gen("bc", TRUE, FALSE, integer(0)), 3:4)

})

test_that("backcross nrec works", {

    # autosome
    expect_equal(test_nrec("bc", 1, 1, FALSE, FALSE, integer(0)), 0)
    expect_equal(test_nrec("bc", 1, 2, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("bc", 2, 1, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("bc", 2, 2, FALSE, FALSE, integer(0)), 0)
    # X female
    expect_equal(test_nrec("bc", 1, 1, TRUE, TRUE, integer(0)), 0)
    expect_equal(test_nrec("bc", 1, 2, TRUE, TRUE, integer(0)), 1)
    expect_equal(test_nrec("bc", 2, 1, TRUE, TRUE, integer(0)), 1)
    expect_equal(test_nrec("bc", 2, 2, TRUE, TRUE, integer(0)), 0)
    # X male
    expect_equal(test_nrec("bc", 3, 3, TRUE, FALSE, integer(0)), 0)
    expect_equal(test_nrec("bc", 3, 4, TRUE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("bc", 4, 3, TRUE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("bc", 4, 4, TRUE, FALSE, integer(0)), 0)

})

test_that("backcross init works", {

    # autosome
    expect_equal(test_init("bc", 1, FALSE, FALSE, integer(0)), log(0.5))
    expect_equal(test_init("bc", 2, FALSE, FALSE, integer(0)), log(0.5))
    # X female
    expect_equal(test_init("bc", 1, TRUE, TRUE, integer(0)), log(0.5))
    expect_equal(test_init("bc", 2, TRUE, TRUE, integer(0)), log(0.5))
    # X male
    expect_equal(test_init("bc", 3, TRUE, FALSE, integer(0)), log(0.5))
    expect_equal(test_init("bc", 4, TRUE, FALSE, integer(0)), log(0.5))

})

test_that("backcross emit works", {

    # autosome
    eps <- 0.01
    expect_equal(test_emit("bc", 0, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 2, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 1, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))
    expect_equal(test_emit("bc", 1, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 2, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 2, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))
    # X female
    expect_equal(test_emit("bc", 0, 1, eps, integer(0), TRUE, TRUE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 2, eps, integer(0), TRUE, TRUE, integer(0)), 0)
    expect_equal(test_emit("bc", 1, 1, eps, integer(0), TRUE, TRUE, integer(0)), log(1-eps))
    expect_equal(test_emit("bc", 1, 2, eps, integer(0), TRUE, TRUE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 2, 1, eps, integer(0), TRUE, TRUE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 2, 2, eps, integer(0), TRUE, TRUE, integer(0)), log(1-eps))
    # X male
    expect_equal(test_emit("bc", 0, 3, eps, integer(0), TRUE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 4, eps, integer(0), TRUE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 1, 3, eps, integer(0), TRUE, FALSE, integer(0)), log(1-eps))
    expect_equal(test_emit("bc", 1, 4, eps, integer(0), TRUE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 3, 3, eps, integer(0), TRUE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("bc", 3, 4, eps, integer(0), TRUE, FALSE, integer(0)), log(1-eps))

    # errors
    expect_equal(test_emit("bc", 0, 0, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 3, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 3, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    # X female
    expect_equal(test_emit("bc", 0, 0, eps, integer(0), TRUE, TRUE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 3, eps, integer(0), TRUE, TRUE, integer(0)), 0)
    expect_equal(test_emit("bc", 3, 1, eps, integer(0), TRUE, TRUE, integer(0)), 0)
    # X male
    expect_equal(test_emit("bc", 0, 0, eps, integer(0), TRUE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 1, eps, integer(0), TRUE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 0, 2, eps, integer(0), TRUE, FALSE, integer(0)), 0)
    expect_equal(test_emit("bc", 2, 1, eps, integer(0), TRUE, FALSE, integer(0)), 0)

})

test_that("backcross step works", {

    # autosome
    rf <- 0.01
    expect_equal(test_step("bc", 1, 1, rf, FALSE, FALSE, integer(0)), log(1-rf))
    expect_equal(test_step("bc", 1, 2, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("bc", 2, 1, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("bc", 2, 2, rf, FALSE, FALSE, integer(0)), log(1-rf))
    # X female
    expect_equal(test_step("bc", 1, 1, rf, TRUE, TRUE, integer(0)), log(1-rf))
    expect_equal(test_step("bc", 1, 2, rf, TRUE, TRUE, integer(0)), log(rf))
    expect_equal(test_step("bc", 2, 1, rf, TRUE, TRUE, integer(0)), log(rf))
    expect_equal(test_step("bc", 2, 2, rf, TRUE, TRUE, integer(0)), log(1-rf))
    # X male
    expect_equal(test_step("bc", 3, 3, rf, TRUE, FALSE, integer(0)), log(1-rf))
    expect_equal(test_step("bc", 3, 4, rf, TRUE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("bc", 4, 3, rf, TRUE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("bc", 4, 4, rf, TRUE, FALSE, integer(0)), log(1-rf))

})

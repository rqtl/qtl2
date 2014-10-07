
context("basic HMM functions in doubled haploids")

test_that("doubled haploids check_geno works", {

    # autosome
    expect_true(test_check_geno("dh", 0, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("dh", 1, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("dh", 2, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("dh", 1, FALSE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("dh", 2, FALSE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("dh", 3, TRUE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("dh", 0, FALSE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("dh", 3, FALSE, FALSE, FALSE, numeric(0)))
})

test_that("doubled haploids n_gen works", {

    expect_equal(test_ngen("dh", FALSE), 2)

})

test_that("doubled haploids possible_gen works", {

    expect_equal(test_possible_gen("dh", FALSE, FALSE, numeric(0)), 1:2)
})

test_that("doubled haploids nrec works", {

    # autosome
    expect_equal(test_nrec("dh", 1, 1, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_nrec("dh", 1, 2, FALSE, FALSE, numeric(0)), 1)
    expect_equal(test_nrec("dh", 2, 1, FALSE, FALSE, numeric(0)), 1)
    expect_equal(test_nrec("dh", 2, 2, FALSE, FALSE, numeric(0)), 0)

    # some errors
    expect_error(test_nrec("dh", 0, 1, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("dh", 1, 0, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("dh", 3, 1, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("dh", 1, 3, FALSE, FALSE, numeric(0)))

})

test_that("doubled haploids init works", {

    expect_equal(test_init("dh", 1, FALSE, FALSE, numeric(0)), log(0.5))
    expect_equal(test_init("dh", 2, FALSE, FALSE, numeric(0)), log(0.5))

    # errors
    expect_error(test_init("dh", 0, FALSE, FALSE, numeric(0)))
    expect_error(test_init("dh", 3, FALSE, FALSE, numeric(0)))

})

test_that("doubled haploids emit works", {

    eps <- 0.01
    expect_equal(test_emit("dh", 0, 1, eps, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_emit("dh", 0, 2, eps, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_emit("dh", 1, 1, eps, FALSE, FALSE, numeric(0)), log(1-eps))
    expect_equal(test_emit("dh", 1, 2, eps, FALSE, FALSE, numeric(0)), log(eps))
    expect_equal(test_emit("dh", 2, 1, eps, FALSE, FALSE, numeric(0)), log(eps))
    expect_equal(test_emit("dh", 2, 2, eps, FALSE, FALSE, numeric(0)), log(1-eps))

    # errors
    expect_error(test_emit("dh", 0, 0, eps, FALSE, FALSE, numeric(0)))
    expect_error(test_emit("dh", 0, 3, eps, FALSE, FALSE, numeric(0)))
    expect_error(test_emit("dh", 3, 1, eps, FALSE, FALSE, numeric(0)))

})

test_that("doubled haploids step works", {

    rf <- 0.01
    expect_equal(test_step("dh", 1, 1, rf, FALSE, FALSE, numeric(0)), log(1-rf))
    expect_equal(test_step("dh", 1, 2, rf, FALSE, FALSE, numeric(0)), log(rf))
    expect_equal(test_step("dh", 2, 1, rf, FALSE, FALSE, numeric(0)), log(rf))
    expect_equal(test_step("dh", 2, 2, rf, FALSE, FALSE, numeric(0)), log(1-rf))

    # errors
    expect_error(test_step("dh", 0, 1, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("dh", 1, 0, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("dh", 3, 1, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("dh", 2, 3, rf, FALSE, FALSE, numeric(0)))

})

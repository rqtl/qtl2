
context("basic HMM functions in riself")

test_that("riself check_geno works", {

    expect_true(test_check_geno("riself", 0, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("riself", 1, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("riself", 2, TRUE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("riself", 1, FALSE, FALSE, FALSE, numeric(0)))
    expect_true(test_check_geno("riself", 2, FALSE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("riself", 3, TRUE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("riself", 0, FALSE, FALSE, FALSE, numeric(0)))
    expect_error(test_check_geno("riself", 3, FALSE, FALSE, FALSE, numeric(0)))

})

test_that("riself n_geno works", {
    expect_equal(test_n_geno("riself", FALSE), 2)
})

test_that("riself geno_index works", {

    expect_equal(test_geno_index("riself", FALSE, FALSE, numeric(0)), 0:1)

})

test_that("riself nrec works", {

    expect_equal(test_nrec("riself", 1, 1, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_nrec("riself", 1, 2, FALSE, FALSE, numeric(0)), 1)
    expect_equal(test_nrec("riself", 2, 1, FALSE, FALSE, numeric(0)), 1)
    expect_equal(test_nrec("riself", 2, 2, FALSE, FALSE, numeric(0)), 0)

    # some errors
    expect_error(test_nrec("riself", 0, 1, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("riself", 1, 0, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("riself", 3, 1, FALSE, FALSE, numeric(0)))
    expect_error(test_nrec("riself", 2, 3, FALSE, FALSE, numeric(0)))

})

test_that("riself init works", {

    expect_equal(test_init("riself", 1, FALSE, FALSE, numeric(0)), log(0.5))
    expect_equal(test_init("riself", 2, FALSE, FALSE, numeric(0)), log(0.5))

    # errors
    expect_error(test_init("riself", 0, FALSE, FALSE, numeric(0)))
    expect_error(test_init("riself", 3, FALSE, FALSE, numeric(0)))

})

test_that("riself emit works", {

    eps <- 0.01
    expect_equal(test_emit("riself", 0, 1, eps, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_emit("riself", 0, 2, eps, FALSE, FALSE, numeric(0)), 0)
    expect_equal(test_emit("riself", 1, 1, eps, FALSE, FALSE, numeric(0)), log(1-eps))
    expect_equal(test_emit("riself", 1, 2, eps, FALSE, FALSE, numeric(0)), log(eps))
    expect_equal(test_emit("riself", 2, 1, eps, FALSE, FALSE, numeric(0)), log(eps))
    expect_equal(test_emit("riself", 2, 2, eps, FALSE, FALSE, numeric(0)), log(1-eps))

    # errors
    expect_error(test_emit("riself", 0, 0, eps, FALSE, FALSE, numeric(0)))
    expect_error(test_emit("riself", 0, 3, eps, FALSE, FALSE, numeric(0)))
    expect_error(test_emit("riself", 3, 1, eps, FALSE, FALSE, numeric(0)))

})

test_that("riself step works", {

    rf <- 0.01
    RF <- 2*rf/(1+2*rf)
    expect_equal(test_step("riself", 1, 1, rf, FALSE, FALSE, numeric(0)), log(1-RF))
    expect_equal(test_step("riself", 1, 2, rf, FALSE, FALSE, numeric(0)), log(RF))
    expect_equal(test_step("riself", 2, 1, rf, FALSE, FALSE, numeric(0)), log(RF))
    expect_equal(test_step("riself", 2, 2, rf, FALSE, FALSE, numeric(0)), log(1-RF))

    # errors
    expect_error(test_step("riself", 0, 1, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("riself", 1, 0, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("riself", 3, 1, rf, FALSE, FALSE, numeric(0)))
    expect_error(test_step("riself", 2, 3, rf, FALSE, FALSE, numeric(0)))

})

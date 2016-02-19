context("basic HMM functions in doubled haploids")

test_that("doubled haploids check_geno works", {

    # autosome
    expect_true(test_check_geno("dh", 0, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("dh", 1, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("dh", 2, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("dh", 1, FALSE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("dh", 2, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("dh", 3, TRUE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("dh", 0, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("dh", 3, FALSE, FALSE, FALSE, integer(0)))
})

test_that("doubled haploids n_gen works", {

    expect_equal(test_ngen("dh", FALSE), 2)

})

test_that("doubled haploids possible_gen works", {

    expect_equal(test_possible_gen("dh", FALSE, FALSE, integer(0)), 1:2)
})

test_that("doubled haploids nrec works", {

    # autosome
    expect_equal(test_nrec("dh", 1, 1, FALSE, FALSE, integer(0)), 0)
    expect_equal(test_nrec("dh", 1, 2, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("dh", 2, 1, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("dh", 2, 2, FALSE, FALSE, integer(0)), 0)

})

test_that("doubled haploids init works", {

    expect_equal(test_init("dh", 1, FALSE, FALSE, integer(0)), log(0.5))
    expect_equal(test_init("dh", 2, FALSE, FALSE, integer(0)), log(0.5))

})

test_that("doubled haploids emit works", {

    eps <- 0.01
    expect_equal(test_emit("dh", 0, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("dh", 0, 2, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("dh", 1, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))
    expect_equal(test_emit("dh", 1, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("dh", 2, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("dh", 2, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))

    # errors
    expect_equal(test_emit("dh", 0, 0, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("dh", 0, 3, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("dh", 3, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)

})

test_that("doubled haploids step works", {

    rf <- 0.01
    expect_equal(test_step("dh", 1, 1, rf, FALSE, FALSE, integer(0)), log(1-rf))
    expect_equal(test_step("dh", 1, 2, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("dh", 2, 1, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("dh", 2, 2, rf, FALSE, FALSE, integer(0)), log(1-rf))

})

test_that("geno_names works", {
    expect_equal(geno_names("dh", c("B", "R"), FALSE), c("BB", "RR"))
    expect_equal(geno_names("dh", c("B", "R"), TRUE), c("BB", "RR"))
})

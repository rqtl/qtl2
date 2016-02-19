context("basic HMM functions in haploids")

test_that("haploids check_geno works", {

    # autosome
    expect_true(test_check_geno("haploid", 0, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("haploid", 1, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("haploid", 2, TRUE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("haploid", 1, FALSE, FALSE, FALSE, integer(0)))
    expect_true(test_check_geno("haploid", 2, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("haploid", 3, TRUE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("haploid", 0, FALSE, FALSE, FALSE, integer(0)))
    expect_false(test_check_geno("haploid", 3, FALSE, FALSE, FALSE, integer(0)))
})

test_that("haploids n_gen works", {

    expect_equal(test_ngen("haploid", FALSE), 2)

})

test_that("haploids possible_gen works", {

    expect_equal(test_possible_gen("haploid", FALSE, FALSE, integer(0)), 1:2)
})

test_that("haploids nrec works", {

    # autosome
    expect_equal(test_nrec("haploid", 1, 1, FALSE, FALSE, integer(0)), 0)
    expect_equal(test_nrec("haploid", 1, 2, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("haploid", 2, 1, FALSE, FALSE, integer(0)), 1)
    expect_equal(test_nrec("haploid", 2, 2, FALSE, FALSE, integer(0)), 0)

})

test_that("haploids init works", {

    expect_equal(test_init("haploid", 1, FALSE, FALSE, integer(0)), log(0.5))
    expect_equal(test_init("haploid", 2, FALSE, FALSE, integer(0)), log(0.5))

})

test_that("haploids emit works", {

    eps <- 0.01
    expect_equal(test_emit("haploid", 0, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("haploid", 0, 2, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("haploid", 1, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))
    expect_equal(test_emit("haploid", 1, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("haploid", 2, 1, eps, integer(0), FALSE, FALSE, integer(0)), log(eps))
    expect_equal(test_emit("haploid", 2, 2, eps, integer(0), FALSE, FALSE, integer(0)), log(1-eps))

    # errors
    expect_equal(test_emit("haploid", 0, 0, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("haploid", 0, 3, eps, integer(0), FALSE, FALSE, integer(0)), 0)
    expect_equal(test_emit("haploid", 3, 1, eps, integer(0), FALSE, FALSE, integer(0)), 0)

})

test_that("haploids step works", {

    rf <- 0.01
    expect_equal(test_step("haploid", 1, 1, rf, FALSE, FALSE, integer(0)), log(1-rf))
    expect_equal(test_step("haploid", 1, 2, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("haploid", 2, 1, rf, FALSE, FALSE, integer(0)), log(rf))
    expect_equal(test_step("haploid", 2, 2, rf, FALSE, FALSE, integer(0)), log(1-rf))

})

test_that("geno_names works", {
    expect_equal(geno_names("haploid", c("B", "R"), FALSE), c("B", "R"))
    expect_equal(geno_names("haploid", c("B", "R"), TRUE), c("B", "R"))
})

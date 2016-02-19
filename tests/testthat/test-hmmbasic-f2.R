context("basic HMM functions in intercross")

test_that("intercross check_geno works", {

    # autosome
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0))
    for(i in 1:3)
        expect_true(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0))
    for(i in c(-1, 6))
        expect_false(test_check_geno("f2", i, TRUE, FALSE, FALSE, 0))
    for(i in c(0, 4))
        expect_false(test_check_geno("f2", i, FALSE, FALSE, FALSE, 0))

    # X chromosome female, forward cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0))
    for(i in 1:2)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0))
    for(i in c(-1, 6))
        expect_false(test_check_geno("f2", i, TRUE, TRUE, TRUE, 0))
    for(i in c(0, 3, 4))
        expect_false(test_check_geno("f2", i, FALSE, TRUE, TRUE, 0))

    # X chromosome female, reverse cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1))
    for(i in 3:4)
        expect_true(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1))
    for(i in c(-1, 6))
        expect_false(test_check_geno("f2", i, TRUE, TRUE, TRUE, 1))
    for(i in c(0, 1, 2, 5, 6))
        expect_false(test_check_geno("f2", i, FALSE, TRUE, TRUE, 1))

    # X chromosome male
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0))
    for(i in c(5,6))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0))
    for(i in c(-1, 6))
        expect_false(test_check_geno("f2", i, TRUE, TRUE, FALSE, 0))
    for(i in c(0:4, 7))
        expect_false(test_check_geno("f2", i, FALSE, TRUE, FALSE, 0))

    # X chromosome male reverse cross
    for(i in 0:5)
        expect_true(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1))
    for(i in c(5,6))
        expect_true(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1))
    for(i in c(-1, 6))
        expect_false(test_check_geno("f2", i, TRUE, TRUE, FALSE, 1))
    for(i in c(0:4, 7))
        expect_false(test_check_geno("f2", i, FALSE, TRUE, FALSE, 1))

})

test_that("intercross n_gen works", {

    expect_equal(test_ngen("f2", FALSE), 3)
    expect_equal(test_ngen("f2", TRUE),  6)

})

test_that("intercross possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("f2", FALSE, FALSE, 0), 1:3)
    # X female forward
    expect_equal(test_possible_gen("f2", TRUE, TRUE, 0), 1:2)
    # X female reverse
    expect_equal(test_possible_gen("f2", TRUE, TRUE, 1), 3:4)
    # X male
    expect_equal(test_possible_gen("f2", TRUE, FALSE, 0), 5:6)
    # X male reverse
    expect_equal(test_possible_gen("f2", TRUE, FALSE, 1), 5:6)

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
    for(i in 3:4)
        expect_equal(test_init("f2", i, TRUE, TRUE, 1), log(0.5))
    # X male
    for(i in 5:6)
        expect_equal(test_init("f2", i, TRUE, FALSE, 0), log(0.5))
    # X male reverse
    for(i in 5:6)
        expect_equal(test_init("f2", i, TRUE, FALSE, 1), log(0.5))

})

test_that("intercross emit works", {

    # autosome
    eps <- 0.01
    for(i in 1:3)
        expect_equal(test_emit("f2", 0, i, eps, integer(0), FALSE, FALSE, 0), 0)

    expect_equal(test_emit("f2", 1, 1, eps, integer(0), FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 3, 1, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 4, 1, eps, integer(0), FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2", 5, 1, eps, integer(0), FALSE, FALSE, 0), log(eps))

    expect_equal(test_emit("f2", 1, 2, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 2, 2, eps, integer(0), FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 4, 2, eps, integer(0), FALSE, FALSE, 0), log(1-eps/2))
    expect_equal(test_emit("f2", 5, 2, eps, integer(0), FALSE, FALSE, 0), log(1-eps/2))

    expect_equal(test_emit("f2", 1, 3, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 2, 3, eps, integer(0), FALSE, FALSE, 0), log(eps/2))
    expect_equal(test_emit("f2", 3, 3, eps, integer(0), FALSE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 3, eps, integer(0), FALSE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 3, eps, integer(0), FALSE, FALSE, 0), log(1-eps/2))

    # X female forward
    for(i in 1:2)
        expect_equal(test_emit("f2", 0, i, eps, integer(0), TRUE, TRUE, 0), 0)

    expect_equal(test_emit("f2", 1, 1, eps, integer(0), TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 1, eps, integer(0), TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 3, 1, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 4, 1, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 1, eps, integer(0), TRUE, TRUE, 0), log(eps))

    expect_equal(test_emit("f2", 1, 2, eps, integer(0), TRUE, TRUE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 2, eps, integer(0), TRUE, TRUE, 0), log(1-eps))
    expect_equal(test_emit("f2", 3, 2, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 4, 2, eps, integer(0), TRUE, TRUE, 0), 0)
    expect_equal(test_emit("f2", 5, 2, eps, integer(0), TRUE, TRUE, 0), log(1-eps))

    # X female reverse
    for(i in 3:4)
        expect_equal(test_emit("f2", 0, i, eps, integer(0), TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 1, 4, eps, integer(0), TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 2, 4, eps, integer(0), TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 3, 4, eps, integer(0), TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 4, 4, eps, integer(0), TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 5, 4, eps, integer(0), TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 1, 3, eps, integer(0), TRUE, TRUE, 1), 0)
    expect_equal(test_emit("f2", 2, 3, eps, integer(0), TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 3, 3, eps, integer(0), TRUE, TRUE, 1), log(eps))
    expect_equal(test_emit("f2", 4, 3, eps, integer(0), TRUE, TRUE, 1), log(1-eps))
    expect_equal(test_emit("f2", 5, 3, eps, integer(0), TRUE, TRUE, 1), 0)

    # X male
    for(i in 5:6)
        expect_equal(test_emit("f2", 0, i, eps, integer(0), TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 1, 5, eps, integer(0), TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 2, 5, eps, integer(0), TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 3, 5, eps, integer(0), TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 4, 5, eps, integer(0), TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 5, 5, eps, integer(0), TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 1, 6, eps, integer(0), TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 2, 6, eps, integer(0), TRUE, FALSE, 0), 0)
    expect_equal(test_emit("f2", 3, 6, eps, integer(0), TRUE, FALSE, 0), log(1-eps))
    expect_equal(test_emit("f2", 4, 6, eps, integer(0), TRUE, FALSE, 0), log(eps))
    expect_equal(test_emit("f2", 5, 6, eps, integer(0), TRUE, FALSE, 0), log(1-eps))

    # X male reverse
    for(i in 5:6)
        expect_equal(test_emit("f2", 0, i, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 1, 5, eps, integer(0), TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 2, 5, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 3, 5, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 4, 5, eps, integer(0), TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 5, 5, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 1, 6, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 2, 6, eps, integer(0), TRUE, FALSE, 1), 0)
    expect_equal(test_emit("f2", 3, 6, eps, integer(0), TRUE, FALSE, 1), log(1-eps))
    expect_equal(test_emit("f2", 4, 6, eps, integer(0), TRUE, FALSE, 1), log(eps))
    expect_equal(test_emit("f2", 5, 6, eps, integer(0), TRUE, FALSE, 1), log(1-eps))

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
    expect_equal(test_step("f2", 3, 4, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2", 4, 3, rf, TRUE, TRUE, 1), log(rf))
    expect_equal(test_step("f2", 4, 4, rf, TRUE, TRUE, 1), log(1-rf))
    # X male
    expect_equal(test_step("f2", 5, 5, rf, TRUE, FALSE, 0), log(1-rf))
    expect_equal(test_step("f2", 5, 6, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 6, 5, rf, TRUE, FALSE, 0), log(rf))
    expect_equal(test_step("f2", 6, 6, rf, TRUE, FALSE, 0), log(1-rf))
    # X male reverse
    expect_equal(test_step("f2", 5, 5, rf, TRUE, FALSE, 1), log(1-rf))
    expect_equal(test_step("f2", 5, 6, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2", 6, 5, rf, TRUE, FALSE, 1), log(rf))
    expect_equal(test_step("f2", 6, 6, rf, TRUE, FALSE, 1), log(1-rf))

})

test_that("geno_names works", {
    expect_equal(geno_names("f2", c("B", "R"), FALSE), c("BB", "BR", "RR"))
    expect_equal(geno_names("f2", c("B", "R"), TRUE), c("BB", "BR", "RR", "BY", "RY"))
})

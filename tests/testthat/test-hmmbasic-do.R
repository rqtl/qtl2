context("basic HMM functions in Diversity Outcross")

test_that("DO check_geno works", {

    # observed genotypes
    for(i in 0:5) {
        # Autosome
        expect_true(test_check_geno("do", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("do", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_true(test_check_geno("do", i, TRUE, TRUE, FALSE, 20))
    }
    for(i in c(-1, 6)) {
        # Autosome
        expect_false(test_check_geno("do", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("do", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_false(test_check_geno("do", i, TRUE, TRUE, FALSE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 1:36) {
        # Autosome
        expect_true(test_check_geno("do", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("do", i, FALSE, TRUE, TRUE, 20))
    }
    for(i in c(0, 37)) {
        # Autosome
        expect_false(test_check_geno("do", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("do", i, FALSE, TRUE, TRUE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 36 + 1:8) {
        expect_true(test_check_geno("do", i, FALSE, TRUE, FALSE, 20))
    }
    for(i in c(0:36, 36+9)) {
        expect_false(test_check_geno("do", i, FALSE, TRUE, FALSE, 20))
    }

})

test_that("DO n_gen works", {

    expect_equal(test_ngen("do", FALSE), 36)
    expect_equal(test_ngen("do", TRUE),  36+8)

})

test_that("DO possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("do", FALSE, FALSE, 20), 1:36)

    # X female
    expect_equal(test_possible_gen("do", TRUE, TRUE, 20), 1:36)

    # X male
    expect_equal(test_possible_gen("do", TRUE, FALSE, 20), 36+(1:8))

})

test_that("DO init works", {

    hom <- cumsum(1:8)
    het <- (1:36)[!((1:36) %in% cumsum(1:8))]
    male <- 36 + (1:8)

    for(i in hom) {
        # autosome and female X
        expect_equal(test_init("do", i, FALSE, FALSE, 20), log(1/64))
        expect_equal(test_init("do", i, TRUE,  TRUE,  20), log(1/64))
    }

    for(i in het) {
        # autosome and female X
        expect_equal(test_init("do", i, FALSE, FALSE, 20), log(1/32))
        expect_equal(test_init("do", i, TRUE,  TRUE,  20), log(1/32))
    }

    for(i in male)
        expect_equal(test_init("do", i, TRUE,  FALSE, 20), log(1/8))

})


test_that("DO emit works", {


})

test_that("DO step works", {


})

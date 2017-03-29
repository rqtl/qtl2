context("compare_geno")

test_that("compare_geno works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[1:5, c(18:19,"X")]

    cg <- compare_geno(iron)
    expected <- structure(c(6, 0, 6, 6, 0, 0, 0, 0, 0, 0, 2, 0, 6, 6, 0, 1, 0,
                            4, 6, 0, 0, 0, 0, 0, 0), .Dim = c(5L, 5L),
                          .Dimnames = list(c("1", "2", "3", "4", "5"), c("1", "2", "3", "4", "5")),
                          class = c("compare_geno", "matrix"))
    expect_equal(cg, expected)

    cg_X <- compare_geno(iron[,"X"])
    expectedX <- structure(c(2, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0,
                             2, 2, 0, 0, 0, 0, 0, 0), .Dim = c(5L, 5L),
                           .Dimnames = list(c("1", "2", "3", "4", "5"), c("1", "2", "3", "4", "5")),
                           class = c("compare_geno", "matrix"))
    expect_equal(cg_X, expectedX)

    cg_noX <- compare_geno(iron, omit_x=TRUE)
    expect_equal(cg_noX, expected - expectedX)

})

test_that("compare_geno works when multi-core", {

    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[1:5, c(18:19,"X")]

    cg <- compare_geno(iron)
    cg_mc <- compare_geno(iron, cores=3)
    expect_equal(cg, cg_mc)

})

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

    # test summary()
    expected <- structure(list(ind1 = character(0), ind2 = character(0), prop_match = numeric(0),
                               n_mismatch = numeric(0), n_typed = numeric(0), n_match = numeric(0)),
                          .Names = c("ind1", "ind2", "prop_match", "n_mismatch", "n_typed", "n_match"),
                          row.names = integer(0), class = c("summary.compare_geno", "data.frame"),
                          threshold = 0.9)
    expect_equal(summary(cg), expected)

    expected <- structure(list(ind1 = "3", ind2 = "4", prop_match = 0.666666666666667,
                               n_mismatch = 2, n_typed = 6, n_match = 4),
                          .Names = c("ind1", "ind2", "prop_match", "n_mismatch", "n_typed", "n_match"),
                          row.names = 1L, class = c("summary.compare_geno", "data.frame"),
                          threshold = 0.5)
    expect_equal(summary(cg, 0.5), expected)

})

test_that("compare_geno works when multi-core", {

    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[1:5, c(18:19,"X")]

    cg <- compare_geno(iron)
    cg_mc <- compare_geno(iron, cores=3)
    expect_equal(cg, cg_mc)

})

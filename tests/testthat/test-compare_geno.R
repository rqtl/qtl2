context("compare_geno")

test_that("compare_geno works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[1:5, c(18:19,"X")]

    cg <- compare_geno(iron)
    expected <- structure(c(6, 0, 6, 6, 0, NA, 0, 0, 0, 0, 0.333333333333333,
                            NA, 6, 6, 0, 0.166666666666667, NA, 0.666666666666667, 6, 0,
                            NA, NA, NA, NA, 0), .Dim = c(5L, 5L),
                          .Dimnames = list(c("1", "2", "3", "4", "5"), c("1", "2", "3", "4", "5")),
                          proportion=TRUE,
                          class = c("compare_geno", "matrix"))
    expect_equal(cg, expected)

    cg_X <- compare_geno(iron[,"X"])
    expectedX <- structure(c(2, 0, 2, 2, 0, NA, 0, 0, 0, 0, 0, NA, 2, 2, 0, 0,
                             NA, 1, 2, 0, NA, NA, NA, NA, 0), .Dim = c(5L, 5L),
                           .Dimnames = list(c("1", "2", "3", "4", "5"), c("1", "2", "3", "4", "5")),
                           proportion=TRUE,
                           class = c("compare_geno", "matrix"))
    expect_equal(cg_X, expectedX)

    cg_noX <- compare_geno(iron, omit_x=TRUE)
    expected_noX <- structure(c(4, 0, 4, 4, 0, NA, 0, 0, 0, 0, 0.5, NA, 4, 4, 0,
                                0.25, NA, 0.5, 4, 0, NA, NA, NA, NA, 0), .Dim = c(5L, 5L),
                              .Dimnames = list(c("1", "2", "3", "4", "5"), c("1", "2", "3", "4", "5")),
                              proportion=TRUE,
                              class = c("compare_geno", "matrix"))
    expect_equal(cg_noX, expected_noX)

    # if proportion=FALSE, total is sum of X and autosomes
    expect_equal(compare_geno(iron, proportion=FALSE),
                 compare_geno(iron[,"X"], proportion=FALSE) +
                 compare_geno(iron, omit_x=TRUE, proportion=FALSE))

    # test summary()
    expected <- structure(list(ind1 = character(0), ind2 = character(0), prop_match = numeric(0),
                               n_mismatch = numeric(0), n_typed = numeric(0), n_match = numeric(0),
                               index1=numeric(0), index2=numeric(0)),
                          .Names = c("ind1", "ind2", "prop_match", "n_mismatch", "n_typed", "n_match", "index1", "index2"),
                          row.names = integer(0), class = c("summary.compare_geno", "data.frame"),
                          threshold = 0.9)
    expect_equal(summary(cg), expected)

    expected <- structure(list(ind1 = "3", ind2 = "4", prop_match = 0.666666666666667,
                               n_mismatch = 2, n_typed = 6, n_match = 4, index1=3, index2=4),
                          .Names = c("ind1", "ind2", "prop_match", "n_mismatch", "n_typed", "n_match", "index1", "index2"),
                          row.names = 1L, class = c("summary.compare_geno", "data.frame"),
                          threshold = 0.5)
    expect_equal(summary(cg, 0.5), expected)

    # test max()
    attr(expected, "threshold") <- NULL
    expect_equal(max(cg), expected)

    # summary and max should give the same results whether you use proportion=TRUE or FALSE
    expect_equal(summary(compare_geno(iron, proportion=TRUE)),
                 summary(compare_geno(iron, proportion=FALSE)))
    expect_equal(summary(compare_geno(iron, proportion=TRUE), threshold=0.7),
                 summary(compare_geno(iron, proportion=FALSE), threshold=0.7))
    expect_equal(max(compare_geno(iron, proportion=TRUE)),
                 max(compare_geno(iron, proportion=FALSE)))
    expect_equal(summary(compare_geno(iron, omit_x=TRUE, proportion=TRUE)),
                 summary(compare_geno(iron, omit_x=TRUE, proportion=FALSE)))
    expect_equal(summary(compare_geno(iron, omit_x=TRUE, proportion=TRUE), threshold=0.7),
                 summary(compare_geno(iron, omit_x=TRUE, proportion=FALSE), threshold=0.7))
    expect_equal(max(compare_geno(iron, omit_x=TRUE, proportion=TRUE)),
                 max(compare_geno(iron, omit_x=TRUE, proportion=FALSE)))

})

test_that("compare_geno works when multi-core", {

    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[1:5, c(18:19,"X")]

    cg <- compare_geno(iron)
    cg_mc <- compare_geno(iron, cores=2)
    expect_equal(cg, cg_mc)

})

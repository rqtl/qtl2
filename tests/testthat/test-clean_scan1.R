context("clean scan1 output")

test_that("clean_scan1 works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c("2", "8", "9")]
    pr <- calc_genoprob(iron)
    out <- scan1(pr, iron$pheno)

    # add some messiness
    out[5,1] <- -out[5,1]
    out[7,2] <- -out[7,2]
    out[14,] <- NA

    # expected values:
    expected <- out[-14,]
    expected[5,1] <- NA
    expected[7,2] <- NA
    class(expected) <- class(out)
    attr(expected, "sample_size") <- attr(out, "sample_size")

    expect_equal(clean_scan1(out), expected)
    expect_equal(clean(out), expected)

    expect_equal(clean_scan1(expected), expected)
    expect_equal(clean(expected), expected)

})

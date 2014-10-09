
context("input/output")

test_that("can read grav2 data", {

    zip_file <- system.file("extdata", "grav2.zip", package="qtl2")
    suppressMessages(grav2 <- read_cross2(zip_file))

    # calculate QTL genotype probabilities
    pr <- calc_genoprob(grav2, step=1)

})

test_that("can read iron data", {

    zip_file <- system.file("extdata", "iron.zip", package="qtl2")
    suppressMessages(iron <- read_cross2(zip_file))

    # calculate QTL genotype probabilities
    pr <- calc_genoprob(iron, step=1)

})

context("matrix to strata")

test_that("mat2strata works", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    Xcovar <- get_x_covar(iron)[c(1:3, 146:148, 215:217),]

    expect_equal(mat2strata(Xcovar),
                 c("1"="1|0", "2"="1|0", "3"="1|0",
                   "146"="0|1", "147"="0|1", "148"="0|1",
                   "215"="0|0", "216"="0|0", "217"="0|0"))

})

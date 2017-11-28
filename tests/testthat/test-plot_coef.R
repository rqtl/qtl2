context("plot_coef")

test_that("plot_coef works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,2]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    pheno <- iron$pheno[,1,drop=FALSE]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)

    coef <- scan1coef(probs, pheno, addcovar=covar)

    test_plot_coef <- function() plot(coef, map, columns=1:3)

    vdiffr::expect_doppelganger("plot_coef", test_plot_coef)

})

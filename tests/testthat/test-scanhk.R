context("genome scan by Haley-Knott")
library(qtl)

test_that("genome scan by Haley-Knott same as R/qtl", {

    # data for chr 1
    data(hyper)
    hyper <- hyper[1,]
    hyper2 <- qtl2geno::convert2cross2(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk")

    # inputs for R/qtl2
    pr <- qtl2geno::calc_genoprob(hyper2, step=1)[[1]][,1,,drop=FALSE]
    y <- hyper2$pheno[,1]

    # center
    intercept <- cbind(rep(1, qtl2geno::n_ind(hyper2)))
    pr <- qtl2scan:::calc_resid_linreg_3d(intercept, pr)
    y <- qtl2scan:::calc_resid_linreg(intercept, as.matrix(y))

    # scan
#    rss <- scan_hk_onechr_nocovar(pr,


})

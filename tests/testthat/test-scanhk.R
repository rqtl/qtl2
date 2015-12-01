context("genome scan by Haley-Knott")
library(qtl)

test_that("genome scan by Haley-Knott same as R/qtl", {

    # data for chr 4
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- qtl2geno::convert2cross2(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk")
    lod0 <- out[,3]

    # inputs for R/qtl2
    pr <- qtl2geno::calc_genoprob(hyper2, step=1)[[1]][,2,,drop=FALSE]
    y <- hyper2$pheno[,1]

    # scan
    rss1 <- scan_hk_onechr(pr, as.matrix(y), as.matrix(rep(1, length(y))))
    lod1 <- length(y)/2 * (log10(sum((y-mean(y))^2)) - log10(rss1))
    lod1 <- as.numeric(lod1)

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # try it by centering first
    intercept <- cbind(rep(1, qtl2geno::n_ind(hyper2)))
    pr <- calc_resid_linreg_3d(intercept, pr)
    y <- calc_resid_linreg(intercept, as.matrix(y))

    # scan
    rss2 <- scan_hk_onechr_nocovar(pr, y)
    lod2 <- length(y)/2 * (log10(sum(y^2)) - log10(rss2))
    lod2 <- as.numeric(lod2)

    # as expected?
    expect_equal(lod0, lod2)
    expect_equal(rss1, rss2)

})

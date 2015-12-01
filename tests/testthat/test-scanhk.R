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
    n <- length(y)

    # scan
    rss1 <- scan_hk_onechr(pr, as.matrix(y), as.matrix(rep(1, n)))
    lod1 <- n/2 * (log10(sum((y-mean(y))^2)) - log10(rss1))
    lod1 <- as.numeric(lod1)

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # try it by centering first
    intercept <- cbind(rep(1, n))
    pr_r <- calc_resid_linreg_3d(intercept, pr)
    y_r <- calc_resid_linreg(intercept, as.matrix(y))

    # scan
    rss2 <- scan_hk_onechr_nocovar(pr_r, y_r)
    lod2 <- n/2 * (log10(sum(y_r^2)) - log10(rss2))
    lod2 <- as.numeric(lod2)

    # as expected?
    expect_equal(lod0, lod2)
    expect_equal(rss1, rss2)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w)
    lodw0 <- outw[,3]

    rssw1 <- scan_hk_onechr_weighted(pr, as.matrix(y), as.matrix(rep(1, n)), sqrt(w))
    lodw1 <- n/2 * (log10(sum(lm(y ~ 1, weights=w)$resid^2*w)) - log10(rssw1))
    lodw1 <- as.numeric(lodw1)

    rssw2 <- scan_hk_onechr(pr*sqrt(w), as.matrix(y)*sqrt(w), as.matrix(rep(1, n))*sqrt(w))
    lodw2 <- n/2 * (log10(sum(lm(y ~ 1, weights=w)$resid^2*w)) - log10(rssw2))
    lodw2 <- as.numeric(lodw2)

    # as expected?
    expect_equal(lodw0, lodw1)
    expect_equal(lodw1, lodw2)

})

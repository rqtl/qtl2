context("genome scan by Haley-Knott")
library(qtl)

test_that("genome scan by Haley-Knott works", {

    # data for chr 6
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

    ###
    # direct calculation with lm()
    rss3 <- apply(pr, 3, function(a) sum(lm(y ~ a)$resid^2))
    names(rss3) <- NULL
    expect_equal(as.numeric(rss1), rss3)

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

    ###
    # direct calculation with lm()
    rssw3 <- apply(pr, 3, function(a) sum(lm(y ~ a, weights=w)$resid^2*w))
    names(rssw3) <- NULL
    expect_equal(as.numeric(rssw1), rssw3)

})

test_that("genome scan by Haley-Knott with multiple phenotypes works", {

    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- qtl2geno::convert2cross2(hyper)
    n_phe <- 200
    hyper$pheno <- cbind(permute_nvector(n_phe, hyper$pheno[,1]),
                         hyper$pheno[,2,drop=FALSE]) # sex column left at the end

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk", pheno.col=1:n_phe)
    lod0 <- t(out[,-(1:2)])
    dimnames(lod0) <- NULL

    # inputs for R/qtl2
    pr <- qtl2geno::calc_genoprob(hyper2, step=1)[[1]][,2,,drop=FALSE]
    y <- as.matrix(hyper$pheno[,1:n_phe])
    n <- nrow(y)

    # scan
    rss1 <- scan_hk_onechr(pr, y, as.matrix(rep(1, n)))
    lod1 <- n/2 * (log10(colSums((y-colMeans(y))^2)) - log10(rss1))

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # try it by centering first
    intercept <- cbind(rep(1, n))
    pr_r <- calc_resid_linreg_3d(intercept, pr)
    y_r <- calc_resid_linreg(intercept, as.matrix(y))

    # scan
    rss2 <- scan_hk_onechr_nocovar(pr_r, y_r)
    lod2 <- n/2 * (log10(colSums(y_r^2)) - log10(rss2))

    # as expected?
    expect_equal(lod0, lod2)
    expect_equal(rss1, rss2)

    ###
    # direct calculation with lm()
    rss3 <- apply(pr, 3, function(a) colSums(lm(y ~ a)$resid^2))
    dimnames(rss3) <- NULL
    expect_equal(rss1, rss3)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, phe=1:n_phe)
    lodw0 <- t(outw[,-(1:2)])
    dimnames(lodw0) <- NULL

    rssw1 <- scan_hk_onechr_weighted(pr, y, as.matrix(rep(1, n)), sqrt(w))
    lodw1 <- n/2 * (log10(colSums(lm(y~1, weights=w)$resid^2*w)) - log10(rssw1))

    rssw2 <- scan_hk_onechr(pr*sqrt(w), y*sqrt(w), as.matrix(rep(1, n))*sqrt(w))
    lodw2 <- n/2 * (log10(colSums(lm(y~1, weights=w)$resid^2*w)) - log10(rssw2))

    # as expected?
    expect_equal(lodw0, lodw1)
    expect_equal(lodw1, lodw2)

    ###
    # direct calculation with lm()
    rssw3 <- apply(pr, 3, function(a) colSums(lm(y ~ a, weights=w)$resid^2*w))
    dimnames(rssw3) <- NULL
    expect_equal(rssw1, rssw3)

})

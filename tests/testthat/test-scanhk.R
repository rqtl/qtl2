context("chromosome scan by basic Haley-Knott functions")
library(qtl)

test_that("chromosome scan by Haley-Knott works", {

    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk")
    lod0 <- out[,3]

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
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

test_that("chromosome scan by Haley-Knott with multiple phenotypes works", {

    skip_if(isnt_karl(), "this test only run locally")

    set.seed(20151201)
    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)
    n_phe <- 200
    hyper$pheno <- cbind(permute_nvector(n_phe, hyper$pheno[,1]),
                         hyper$pheno[,2,drop=FALSE]) # sex column left at the end

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk", pheno.col=1:n_phe)
    lod0 <- t(out[,-(1:2)])
    dimnames(lod0) <- NULL

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
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


test_that("chromosome scan by Haley-Knott works with additive covariates", {

    skip_if(isnt_karl(), "this test only run locally")

    set.seed(20151201)

    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    p <- hyper$pheno[,1]; p <- (p-min(p))/max(p)
    x <- rbinom(nind(hyper), 1, prob=p)
    out <- scanone(hyper, addcovar=x, method="hk")
    lod0 <- out[,3]

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
    y <- hyper2$pheno[,1]
    n <- length(y)

    # scan
    rss1 <- scan_hk_onechr(pr, as.matrix(y), cbind(1, x))
    lod1 <- n/2 * (log10(sum(lm(y~x)$resid^2)) - log10(rss1))
    lod1 <- as.numeric(lod1)

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # direct calculation with lm()
    rss2 <- apply(pr, 3, function(a) sum(lm(y ~ x+a)$resid^2))
    names(rss2) <- NULL
    expect_equal(as.numeric(rss1), rss2)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, addcovar=x)
    lodw0 <- outw[,3]

    rssw1 <- scan_hk_onechr_weighted(pr, as.matrix(y), cbind(1,x), sqrt(w))
    lodw1 <- n/2 * (log10(sum(lm(y ~ x, weights=w)$resid^2*w)) - log10(rssw1))
    lodw1 <- as.numeric(lodw1)

    # as expected?
    expect_equal(lodw0, lodw1)

    ###
    # direct calculation with lm()
    rssw2 <- apply(pr, 3, function(a) sum(lm(y ~ x+a, weights=w)$resid^2*w))
    names(rssw2) <- NULL
    expect_equal(as.numeric(rssw1), rssw2)

})

test_that("chromosome scan by Haley-Knott with multiple phenotypes and an additive covariate works", {

    skip_if(isnt_karl(), "this test only run locally")

    set.seed(20151201)
    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)
    n_phe <- 200
    hyper$pheno <- cbind(permute_nvector(n_phe, hyper$pheno[,1]),
                         hyper$pheno[,2,drop=FALSE]) # sex column left at the end

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    p <- hyper$pheno[,1]; p <- (p-min(p))/max(p)
    x <- rbinom(nind(hyper), 1, prob=p)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    lod0 <- t(out[,-(1:2)])
    dimnames(lod0) <- NULL

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
    y <- as.matrix(hyper$pheno[,1:n_phe])
    n <- nrow(y)

    # scan
    rss1 <- scan_hk_onechr(pr, y, cbind(1,x))
    lod1 <- n/2 * (log10(colSums(lm(y~x)$resid^2)) - log10(rss1))

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # direct calculation with lm()
    rss2 <- apply(pr, 3, function(a) colSums(lm(y ~ x+a)$resid^2))
    dimnames(rss2) <- NULL
    expect_equal(rss1, rss2)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, addcovar=x, phe=1:n_phe)
    lodw0 <- t(outw[,-(1:2)])
    dimnames(lodw0) <- NULL

    rssw1 <- scan_hk_onechr_weighted(pr, y, cbind(1,x), sqrt(w))
    lodw1 <- n/2 * (log10(colSums(lm(y~x, weights=w)$resid^2*w)) - log10(rssw1))

    # as expected?
    expect_equal(lodw0, lodw1)

    ###
    # direct calculation with lm()
    rssw2 <- apply(pr, 3, function(a) colSums(lm(y ~ x+a, weights=w)$resid^2*w))
    dimnames(rssw2) <- NULL
    expect_equal(rssw1, rssw2)

})


test_that("chromosome scan by Haley-Knott works with interactive covariates", {

    skip_if(isnt_karl(), "this test only run locally")

    set.seed(20151201)

    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    p <- hyper$pheno[,1]; p <- (p-min(p))/max(p)
    x <- rbinom(nind(hyper), 1, prob=p)
    out <- scanone(hyper, addcovar=x, intcovar=x, method="hk")
    lod0 <- out[,3]

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
    y <- hyper2$pheno[,1]
    n <- length(y)

    # scan
    rss1 <- scan_hk_onechr_intcovar_highmem(pr, as.matrix(y), cbind(1, x), as.matrix(x))
    rss1_lm <- scan_hk_onechr_intcovar_lowmem(pr, as.matrix(y), cbind(1, x), as.matrix(x))
    expect_equal(rss1, rss1_lm)
    lod1 <- n/2 * (log10(sum(lm(y~x)$resid^2)) - log10(rss1))
    lod1 <- as.numeric(lod1)

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # direct calculation with lm()
    rss2 <- apply(pr, 3, function(a) sum(lm(y ~ x*a)$resid^2))
    names(rss2) <- NULL
    expect_equal(as.numeric(rss1), rss2)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, addcovar=x, intcovar=x)
    lodw0 <- outw[,3]

    rssw1 <- scan_hk_onechr_intcovar_weighted_highmem(pr, as.matrix(y), cbind(1,x),
                                                      as.matrix(x), sqrt(w))
    rssw1_lm <- scan_hk_onechr_intcovar_weighted_lowmem(pr, as.matrix(y), cbind(1,x),
                                                        as.matrix(x), sqrt(w))
    expect_equal(rssw1, rssw1_lm)
    lodw1 <- n/2 * (log10(sum(lm(y ~ x, weights=w)$resid^2*w)) - log10(rssw1))
    lodw1 <- as.numeric(lodw1)

    # as expected?
    expect_equal(lodw0, lodw1)

    ###
    # direct calculation with lm()
    rssw2 <- apply(pr, 3, function(a) sum(lm(y ~ x*a, weights=w)$resid^2*w))
    names(rssw2) <- NULL
    expect_equal(as.numeric(rssw1), rssw2)

})

test_that("chromosome scan by Haley-Knott with multiple phenotypes and an interactive covariate works", {

    skip_if(isnt_karl(), "this test only run locally")

    set.seed(20151201)
    # data for chr 6
    data(hyper)
    hyper <- hyper[6,]
    hyper2 <- convert2cross2(hyper)
    n_phe <- 200
    hyper$pheno <- cbind(permute_nvector(n_phe, hyper$pheno[,1]),
                         hyper$pheno[,2,drop=FALSE]) # sex column left at the end

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    p <- hyper$pheno[,1]; p <- (p-min(p))/max(p)
    x <- rbinom(nind(hyper), 1, prob=p)
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, pheno.col=1:n_phe)
    lod0 <- t(out[,-(1:2)])
    dimnames(lod0) <- NULL

    # inputs for R/qtl2
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- calc_genoprob(hyper2, map)[[1]][,2,,drop=FALSE]
    y <- as.matrix(hyper$pheno[,1:n_phe])
    n <- nrow(y)

    # scan
    rss1 <- scan_hk_onechr_intcovar_highmem(pr, y, cbind(1,x), as.matrix(x))
    rss1_lm <- scan_hk_onechr_intcovar_lowmem(pr, y, cbind(1,x), as.matrix(x))
    expect_equal(rss1, rss1_lm)
    lod1 <- n/2 * (log10(colSums(lm(y~x)$resid^2)) - log10(rss1))

    # as expected?
    expect_equal(lod0, lod1)

    ###
    # direct calculation with lm()
    rss2 <- apply(pr, 3, function(a) colSums(lm(y ~ x*a)$resid^2))
    dimnames(rss2) <- NULL
    expect_equal(rss1, rss2)

    ##############################
    # weighted scan
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, addcovar=x, intcovar=x, phe=1:n_phe)
    lodw0 <- t(outw[,-(1:2)])
    dimnames(lodw0) <- NULL

    rssw1 <- scan_hk_onechr_intcovar_weighted_highmem(pr, y, cbind(1,x), as.matrix(x), sqrt(w))
    rssw1_lm <- scan_hk_onechr_intcovar_weighted_lowmem(pr, y, cbind(1,x), as.matrix(x), sqrt(w))
    expect_equal(rssw1, rssw1_lm)
    lodw1 <- n/2 * (log10(colSums(lm(y~x, weights=w)$resid^2*w)) - log10(rssw1))

    # as expected?
    expect_equal(lodw0, lodw1)

    ###
    # direct calculation with lm()
    rssw2 <- apply(pr, 3, function(a) colSums(lm(y ~ x*a, weights=w)$resid^2*w))
    dimnames(rssw2) <- NULL
    expect_equal(rssw1, rssw2)

})

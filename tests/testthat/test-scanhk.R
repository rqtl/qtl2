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

    # inputs for R/qtl2
    pr <- qtl2geno::calc_genoprob(hyper2, step=1)[[1]][,2,,drop=FALSE]
    y <- hyper2$pheno[,1]

    # center
    intercept <- cbind(rep(1, qtl2geno::n_ind(hyper2)))
    pr <- qtl2scan:::calc_resid_linreg_3d(intercept, pr)
    y <- qtl2scan:::calc_resid_linreg(intercept, as.matrix(y))

    # scan
    rss <- qtl2scan:::scan_hk_onechr_nocovar(pr, y)
    lod <- qtl2geno::n_ind(hyper2)/2 * (log10(sum(y^2)) - log10(rss))
    lod <- as.numeric(lod)

    # as expected?
    expect_equal(out[,3], lod)

})

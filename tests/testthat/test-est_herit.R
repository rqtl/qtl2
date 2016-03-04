context("estimate heritability with LMM by est_herit")


test_that("est_herit works with intercross", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)

    # scan just chr 19; compare est_herit to those
    out_reml <- scan1_lmm(probs[,"19"], iron$pheno, kinship, reml=TRUE)
    out_ml <- scan1_lmm(probs[,"19"], iron$pheno, kinship, reml=FALSE)

    expect_equal(est_herit(iron$pheno, kinship, reml=TRUE)[1:2], # [1:2] to strip off attributes
                 attr(out_reml, "hsq")[1,]) # [1,] to convert to vector

    expect_equal(est_herit(iron$pheno, kinship, reml=FALSE)[1:2], # [1:2] to strip off attributes
                 attr(out_ml, "hsq")[1,]) # [1,] to convert to vector

})

test_that("est_herit with intercross with an additive covariate", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1_lmm(probs[,"19"], iron$pheno, kinship, addcovar=X, reml=TRUE)
    out_ml <- scan1_lmm(probs[,"19"], iron$pheno, kinship, addcovar=X, reml=FALSE)

    expect_equal(est_herit(iron$pheno, kinship, addcovar=X, reml=TRUE)[1:2], # [1:2] to strip off attributes
                 attr(out_reml, "hsq")[1,]) # [1,] to convert to vector

    expect_equal(est_herit(iron$pheno, kinship, addcovar=X, reml=FALSE)[1:2], # [1:2] to strip off attributes
                 attr(out_ml, "hsq")[1,]) # [1,] to convert to vector

})

context("estimate heritability with LMM by est_herit")


test_that("est_herit works with intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)

    # scan just chr 19; compare est_herit to those
    out_reml <- scan1(probs[,"19"], iron$pheno, kinship, reml=TRUE)
    out_ml <- scan1(probs[,"19"], iron$pheno, kinship, reml=FALSE)

    expect_equal(est_herit(iron$pheno, kinship, reml=TRUE)[1:2], # [1:2] to strip off attributes
                 attr(out_reml, "hsq")[1,]) # [1,] to convert to vector

    expect_equal(est_herit(iron$pheno, kinship, reml=FALSE)[1:2], # [1:2] to strip off attributes
                 attr(out_ml, "hsq")[1,]) # [1,] to convert to vector

    # try when some individuals missing from one or the other dataset
    subind <- 1:100
    expected <- est_herit(iron$pheno[subind,,drop=FALSE], kinship[subind,subind])
    expect_equal(est_herit(iron$pheno, kinship[subind, subind]), expected)
    expect_equal(est_herit(iron$pheno[subind,,drop=FALSE], kinship), expected)

})

test_that("est_herit with intercross with an additive covariate", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs[,"19"], iron$pheno, kinship, addcovar=X, reml=TRUE)
    out_ml <- scan1(probs[,"19"], iron$pheno, kinship, addcovar=X, reml=FALSE)

    expect_equal(est_herit(iron$pheno, kinship, addcovar=X, reml=TRUE)[1:2], # [1:2] to strip off attributes
                 attr(out_reml, "hsq")[1,]) # [1,] to convert to vector

    expect_equal(est_herit(iron$pheno, kinship, addcovar=X, reml=FALSE)[1:2], # [1:2] to strip off attributes
                 attr(out_ml, "hsq")[1,]) # [1,] to convert to vector

    # try when some individuals missing from one or the other dataset
    subind <- c(1:50,101:150)
    expected <- est_herit(iron$pheno[subind,,drop=FALSE], kinship[subind,subind], addcovar=X[subind])
    expect_equal(est_herit(iron$pheno, kinship[subind, subind], addcovar=X), expected)
    expect_equal(est_herit(iron$pheno[subind,,drop=FALSE], kinship, addcovar=X), expected)
    expect_equal(est_herit(iron$pheno, kinship, addcovar=X[subind]), expected)
    expect_equal(est_herit(iron$pheno[subind,,drop=FALSE], kinship[subind, subind], addcovar=X), expected)
    expect_equal(est_herit(iron$pheno[subind,,drop=FALSE], kinship, addcovar=X[subind]), expected)
    expect_equal(est_herit(iron$pheno, kinship[subind,subind], addcovar=X[subind]), expected)

})

test_that("est_herit handles dependent covariate columns", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)

    covar2 <- cbind(a=covar,b=0)
    pheno[1:2,1] <- NA
    covar2[1:2,2] <- 1
    hsq <- est_herit(pheno, kinship, covar)
    hsq2 <- est_herit(pheno, kinship, covar2)
    hsq3 <- est_herit(pheno, kinship, covar2[,-2, drop=FALSE])

    expect_equal(hsq[1], hsq2[1])
    expect_equal(hsq2[2], structure(0.545220193658054, .Names = "spleen"), tolerance=5e-7)
    expect_equal(hsq, hsq3)

    pheno[1:2,1:2] <- NA
    hsq <- est_herit(pheno, kinship, covar)
    hsq2 <- est_herit(pheno, kinship, covar2)
    hsq3 <- est_herit(pheno, kinship, covar2[,-2, drop=FALSE])

    expect_equal(hsq, hsq2)
    expect_equal(hsq, hsq3)

})

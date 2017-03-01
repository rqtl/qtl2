context("subset_scan1")

test_that("subset_scan1 works for intercross with two phenotypes", {

    # load qtl2geno package for data and genoprob calculation
    library(qtl2geno)

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # calculate genotype probabilities
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    # grab phenotypes and covariates; ensure that covariates have names attribute
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # perform genome scan
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # subset one column
    expected <- unclass(out)[,2,drop=FALSE]
    attr(expected, "sample_size") <- attr(out, "sample_size")[2]
    class(expected) <- c("scan1", "matrix")
    expect_equal(subset(out, map, lodcolumn=2), expected)
    expect_equal(out[map,,2], expected)

    # subset chr 2, 8, and 9
    keep <- (map2chr(map) %in% c("2", "8", "9"))
    expected <- unclass(out)[keep,,drop=FALSE]
    attr(expected, "sample_size") <- attr(out, "sample_size")
    class(expected) <- c("scan1", "matrix")
    expect_equal(subset(out, map, chr=c("2", "8", "9")), expected)
    expect_equal(out[map, c("2", "8", "9"),], expected)

    # subset chr 2, 8, and 9 and lodcolumn 2
    expected <- unclass(expected)[,2,drop=FALSE]
    attr(expected, "sample_size") <- attr(out, "sample_size")[2]
    class(expected) <- c("scan1", "matrix")
    expect_equal(subset(out, map, c("2", "8", "9"), 2), expected)
    expect_equal(out[map, c("2", "8", "9"),2], expected)

})

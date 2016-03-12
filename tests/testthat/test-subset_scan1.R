context("subset_scan1")

test_that("subset_scan1 works for intercross with two phenotypes", {

    # load qtl2geno package for data and genoprob calculation
    library(qtl2geno)

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # calculate genotype probabilities
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)

    # grab phenotypes and covariates; ensure that covariates have names attribute
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # perform genome scan
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # subset one column
    expected <- out
    expected$lod <- out$lod[,2,drop=FALSE]
    expected$sample_size <- out$sample_size[2]
    expect_equal(subset(out, lodcolumn=2), expected)
    expect_equal(out[,2], expected)

    # subset chr 2, 8, and 9
    index <- split(1:nrow(out$lod), rep(seq(along=out$map), sapply(out$map, length)))
    keep <- unlist(index[c(2, 8, 9)])
    expected <- out
    expected$lod <- out$lod[keep,,drop=FALSE]
    expected$map <- out$map[c("2", "8", "9")]
    expect_equal(subset(out, chr=c("2", "8", "9")), expected)
    expect_equal(out[c("2", "8", "9"),], expected)

    # subset chr 2, 8, and 9 and lodcolumn 2
    expected$lod <- expected$lod[,2,drop=FALSE]
    expected$sample_size <- expected$sample_size[2]
    expect_equal(subset(out, c("2", "8", "9"), 2), expected)
    expect_equal(out[c("2", "8", "9"),2], expected)

})

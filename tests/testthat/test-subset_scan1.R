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
    out <- scan1(probs, pheno, covar, Xcovar)

    attrib <- attributes(out)
    expected <- out[,2,drop=FALSE]
    attrib[["sample_size"]] <- attrib[["sample_size"]][2]
    for(i in seq(along=attrib)) {
        nam <- names(attrib)[i]
        if(nam %in% c("dim", "dimnames")) next
        attr(expected, nam) <- attrib[[i]]
    }
    expect_equal(subset(out, column=2), expected)

})

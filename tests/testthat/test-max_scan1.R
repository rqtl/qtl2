context("max_scan1")

test_that("max_scan1 works for intercross with two phenotypes", {

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

    # maximum of first column
    expected <- data.frame(chr="16",
                           pos=28.6,
                           liver=max(unclass(out[map,,1])),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c16.loc29"
    expect_equal(max(out,map), expected)

    # maximum of spleen column
    expected <- data.frame(chr="9",
                           pos=56.6,
                           spleen=max(unclass(out[map,,2])),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c9.loc57"
    expect_equal(max(out, map, lodcolumn="spleen"), expected)

    # maximum of first column on chr 2
    expected <- data.frame(chr="2",
                           pos=56.8,
                           liver=max(unclass(out[map,"2",1])),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "D2Mit17"
    expect_equal(max(out, map, chr="2"), expected)

})

test_that("maxlod works for intercross with two phenotypes", {

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

    # overall max
    expect_equal(maxlod(out), max(unclass(out)))

    expect_equal(maxlod(out, map, c("2", "9")),
                 max(unclass(out[map, c("2", "9"), ])))

    expect_equal(maxlod(out, map, "2"), max(unclass(out[map,"2",]), na.rm=TRUE))
    expect_equal(maxlod(out, map, "9"), max(unclass(out[map,"9",]), na.rm=TRUE))

    expect_equal( maxlod(out, map, c("2","8", "11")),
                 max( maxlod(out, map, "2"), maxlod(out, map, "8"), maxlod(out, map, "11") ) )

    # really could be using integers
    expect_equal( maxlod(out, map, c(2,8, 11)),
                 max( maxlod(out, map, 2), maxlod(out, map, 8), maxlod(out, map, 11) ) )

})

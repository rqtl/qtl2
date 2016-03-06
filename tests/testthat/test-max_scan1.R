context("max_scan1")

test_that("max_scan1 works for intercross with two phenotypes", {

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

    # maximum of first column
    expected <- data.frame(chr="16",
                           pos=28.6,
                           liver=max(out[,1]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c16.loc29"
    expect_equal(max(out), expected)

    # maximum of spleen column
    expected <- data.frame(chr="9",
                           pos=56.6,
                           spleen=max(out[,2]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c9.loc57"
    expect_equal(max(out, column="spleen"), expected)

    # maximum of first column on chr 2
    expected <- data.frame(chr="2",
                           pos=56.8,
                           liver=max(out[chr_scan1(out)=="2",1]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "D2Mit17"
    expect_equal(max(out, chr="2"), expected)

})

test_that("maxlod works for intercross with two phenotypes", {

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

    # overall max
    expect_equal(maxlod(out), max(unclass(out)))

    expect_equal(maxlod(out, c("2", "9")),
                 max(unclass(subset(out, chr=c("2", "9")))))

    expect_equal(maxlod(out, "2"), max(out[chr_scan1(out)=="2",], na.rm=TRUE))
    expect_equal(maxlod(out, "9"), max(out[chr_scan1(out)=="9",], na.rm=TRUE))

    expect_equal( maxlod(out, c("2","8", "11")),
                 max( maxlod(out, "2"), maxlod(out, "8"), maxlod(out, "11") ) )

    # really could be using integers
    expect_equal( maxlod(out, c(2,8, 11)),
                 max( maxlod(out, 2), maxlod(out, 8), maxlod(out, 11) ) )

})

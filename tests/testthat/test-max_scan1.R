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
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # maximum of first column
    expected <- data.frame(chr="16",
                           pos=28.6,
                           liver=max(out$lod[,1]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c16.loc29"
    expect_equal(max(out), expected)

    # maximum of spleen column
    expected <- data.frame(chr="9",
                           pos=56.6,
                           spleen=max(out$lod[,2]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c9.loc57"
    expect_equal(max(out, lodcolumn="spleen"), expected)

    # maximum of first column on chr 2
    index <- split(1:nrow(out$lod), rep(seq(along=out$map), sapply(out$map, length)))
    keep <- unlist(index[[2]])
    expected <- data.frame(chr="2",
                           pos=56.8,
                           liver=max(out$lod[keep,1]),
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
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # overall max
    expect_equal(maxlod(out), max(out$lod))

    expect_equal(maxlod(out, c("2", "9")),
                 max(subset(out, chr=c("2", "9"))$lod))

    expect_equal(maxlod(out, "2"), max(out["2",]$lod, na.rm=TRUE))
    expect_equal(maxlod(out, "9"), max(out["9",]$lod, na.rm=TRUE))

    expect_equal( maxlod(out, c("2","8", "11")),
                 max( maxlod(out, "2"), maxlod(out, "8"), maxlod(out, "11") ) )

    # really could be using integers
    expect_equal( maxlod(out, c(2,8, 11)),
                 max( maxlod(out, 2), maxlod(out, 8), maxlod(out, 11) ) )

})

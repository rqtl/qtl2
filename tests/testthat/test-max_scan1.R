context("max_scan1")

test_that("max_scan1 works for intercross with two phenotypes", {

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

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
                           liver=max(out[,1]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c16.loc29"
    expect_equal(max(out,map), expected)

    # maximum of spleen column
    expected <- data.frame(chr="9",
                           pos=56.6,
                           spleen=max(out[,2]),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "c9.loc57"
    expect_equal(max(out, map, lodcolumn="spleen"), expected)

    # maximum of first column on chr 2
    expected <- data.frame(chr="2",
                           pos=56.8,
                           liver=maxlod(subset(out, map, "2", 1)),
                           stringsAsFactors=FALSE)
    rownames(expected) <- "D2Mit17"
    expect_equal(max(out, map, chr="2"), expected)

    # test that it works if output or map are subsetted
    expect_equal( max(out, map, chr="2"), max(subset(out, map, chr="2"), map) )
    expect_equal( max(out, map, chr="2"), max(out, map["2"]) )

    # test that it works if output has rows shuffled
    out_shuffled <- out[sample(1:nrow(out)),,drop=FALSE]
    class(out_shuffled) <- c("scan1", "matrix")
    expect_equal( max(out_shuffled, map), max(out, map))

})

test_that("maxlod works for intercross with two phenotypes", {

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

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
                 max(unclass(subset(out, map, c("2", "9")))))

    expect_equal(maxlod(out, map, "2"), max(unclass(subset(out, map,"2")), na.rm=TRUE))
    expect_equal(maxlod(out, map, "9"), max(unclass(subset(out, map,"9")), na.rm=TRUE))

    expect_equal( maxlod(out, map, c("2","8", "11")),
                 max( maxlod(out, map, "2"), maxlod(out, map, "8"), maxlod(out, map, "11") ) )

    # really could be using integers
    expect_equal( maxlod(out, map, c(2,8, 11)),
                 max( maxlod(out, map, 2), maxlod(out, map, 8), maxlod(out, map, 11) ) )

    # test that it works if output or map are subsetted
    expect_equal( max(out, map, chr="2"), max(subset(out, map, chr="2"), map) )
    expect_equal( max(out, map, chr="2"), max(out, map["2"]) )

    # test that it works if output has rows shuffled
    out_shuffled <- out[sample(1:nrow(out)),,drop=FALSE]
    class(out_shuffled) <- c("scan1", "matrix")
    expect_equal( max(out_shuffled, map), max(out, map))

})

context("scan1max")

test_that("scan1max works", {

    iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
    iron <- iron[,c("1", "19", "X")]
    Xc <- get_x_covar(iron)
    X <- Xc[,1]
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    pr <- calc_genoprob(iron, map, error_prob=0.002)
    k <- calc_kinship(pr)
    k_loco <- calc_kinship(pr, "loco")
    chr <- factor(rep(names(map), sapply(map, length)), names(map))

    out <- scan1(pr, iron$pheno)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno), expected)
    expect_equal(scan1max(pr, iron$pheno, by_chr=TRUE), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, by_chr=TRUE), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, addcovar=X)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, by_chr=TRUE), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, kinship=k)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    for(at in c("sample_size", "hsq")) {
        attr(expected, at) <- attr(expected_bychr, at) <- attr(out, at)
    }
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, kinship=k), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, kinship=k, by_chr=TRUE), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    for(at in c("sample_size", "hsq")) {
        attr(expected, at) <- attr(expected_bychr, at) <- attr(out, at)
    }
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco, by_chr=TRUE), expected_bychr)

})


test_that("scan1max works with multicore", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
    iron <- iron[,c("1", "19", "X")]
    Xc <- get_x_covar(iron)
    X <- Xc[,1]
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    pr <- calc_genoprob(iron, map, error_prob=0.002)
    k <- calc_kinship(pr)
    k_loco <- calc_kinship(pr, "loco")
    chr <- factor(rep(names(map), sapply(map, length)), names(map))

    out <- scan1(pr, iron$pheno, cores=4)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno, cores=4), expected)
    expect_equal(scan1max(pr, iron$pheno, by_chr=TRUE, cores=4), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, cores=4)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, cores=4), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, by_chr=TRUE, cores=4), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, addcovar=X, cores=4)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    attr(expected, "sample_size") <- attr(expected_bychr, "sample_size") <- attr(out, "sample_size")
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, cores=4), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, by_chr=TRUE, cores=4), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, kinship=k, cores=4)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    for(at in c("sample_size", "hsq")) {
        attr(expected, at) <- attr(expected_bychr, at) <- attr(out, at)
    }
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, kinship=k, cores=4), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, kinship=k, by_chr=TRUE, cores=4), expected_bychr)

    out <- scan1(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco, cores=4)
    expected <- apply(out, 2, max, na.rm=TRUE)
    expected_bychr <- apply(out, 2, tapply, chr, max, na.rm=TRUE)
    for(at in c("sample_size", "hsq")) {
        attr(expected, at) <- attr(expected_bychr, at) <- attr(out, at)
    }
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco, cores=4), expected)
    expect_equal(scan1max(pr, iron$pheno, Xcovar=Xc, addcovar=X, kinship=k_loco, by_chr=TRUE, cores=4),
                 expected_bychr)

})

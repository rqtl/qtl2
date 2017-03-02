context("BLUP estimates of QTL effects")

calc_blup <-
    function(probs, pheno, kinship=NULL, addcovar=NULL, tol=1e-12,
             quiet=TRUE)
{
    addcovar <- cbind(rep(1, length(pheno)), addcovar)
    if(!is.null(kinship)) {
        Ke <- decomp_kinship(kinship)
        Keval <- Ke$values
        Kevec <- Ke$vectors

        pheno <- Kevec %*% pheno
        addcovar <- Kevec %*% addcovar
        probs <- Kevec %*% probs

        lmmfit <- Rcpp_fitLMM(Keval, pheno, addcovar, tol=tol)
        hsq <- lmmfit$hsq
        if(!quiet) message("hsq_pg: ", hsq)
        wts <- 1/sqrt(hsq * Keval + (1-hsq))

        pheno <- wts * pheno
        addcovar <- wts * addcovar
        probs <- wts * probs
    }

    k <- probs %*% t(probs)
    ke <- decomp_kinship(k)
    keval <- ke$values
    kevec <- ke$vectors

    pheno <- kevec %*% pheno
    addcovar <- kevec %*% addcovar

    lmmfit <- Rcpp_fitLMM(keval, pheno, addcovar, tol=tol)
    beta <- lmmfit$beta
    hsq <- lmmfit$hsq
    if(!quiet) message("hsq_qtl: ", hsq)
    wts <- hsq * keval + (1-hsq)
    u <- as.numeric( t(kevec %*% probs) %*% diag(hsq/wts) %*%  (pheno - addcovar %*% beta) )
    c(u, beta)
}


library(qtl2geno)
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
phe <- iron$pheno[,1,drop=FALSE]

test_that("scan1blup works with no kinship matrix", {

    pr <- calc_genoprob(iron[,"16"])

    blup <- scan1blup(pr, phe)
    blup_se <- scan1blup(pr, phe, se=TRUE)
    expect_equivalent(blup, blup_se)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, preserve_intercept=TRUE)
    expect_equivalent(unclass(blup), unclass(blup_int)[,1:3] + unclass(blup_int)[,4])

    # brute force BLUPs
    for(i in 1:dim(pr[[1]])[[3]]) {
        blup_alt <- calc_blup(pr[[1]][,,i], phe)
        names(blup_alt) <- colnames(blup_int)
        expect_equal(unclass(blup_int)[i,], blup_alt, tol=1e-6)
    }

    # repeat with an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, addcovar=sex)
    blup_se <- scan1blup(pr, phe, addcovar=sex, se=TRUE)
    expect_equivalent(blup, blup_se)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, addcovar=sex, preserve_intercept=TRUE)
    expect_equivalent(unclass(blup), unclass(blup_int)[,c(1:3,5)] + cbind(unclass(blup_int)[,c(4,4,4)], 0))

    # brute force BLUPs
    for(i in 1:dim(pr[[1]])[[3]]) {
        blup_alt <- calc_blup(pr[[1]][,,i], phe, addcovar=sex)
        names(blup_alt) <- colnames(blup_int)
        expect_equal(unclass(blup_int)[i,], blup_alt, tolerance=1e-5)
    }

})

test_that("scan1blup works with kinship matrix", {

    pr <- calc_genoprob(iron)
    K <- calc_kinship(pr[,c(1:15,17:19,"X")])
    pr <- pr[,"16"]

    # sex as an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, K, sex)
    blup_se <- scan1blup(pr, phe, K, sex, se=TRUE)
    expect_equivalent(blup, blup_se)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, K, sex, preserve_intercept=TRUE)
    expect_equivalent(unclass(blup), unclass(blup_int)[,c(1:3,5)] + cbind(unclass(blup_int)[,c(4,4,4)], 0))

    # brute force BLUPs
    for(i in 1:dim(pr[[1]])[[3]]) {
        blup_alt <- calc_blup(pr[[1]][,,i], phe, K, sex)
        names(blup_alt) <- colnames(blup_int)
        expect_equal(unclass(blup_int)[i,], blup_alt, tolerance=1e-5)
    }

})

test_that("scan1blup works with kinship matrix on another chromosome", {

    pr <- calc_genoprob(iron)
    K <- calc_kinship(pr[,c(1:10,12:19,"X")])
    pr <- pr[,"11"]

    # sex as an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, K, sex)
    blup_se <- scan1blup(pr, phe, K, sex, se=TRUE)
    expect_equivalent(blup, blup_se)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, K, sex, preserve_intercept=TRUE)
    expect_equivalent(unclass(blup), unclass(blup_int)[,c(1:3,5)] + cbind(unclass(blup_int)[,c(4,4,4)], 0))

    # brute force BLUPs
    for(i in 1:dim(pr[[1]])[[3]]) {
        blup_alt <- calc_blup(pr[[1]][,,i], phe, K, sex)
        names(blup_alt) <- colnames(blup_int)
        expect_equal(unclass(blup_int)[i,], blup_alt, tolerance=1e-5)
    }

})

test_that("scan1blup deals with mismatching individuals", {
    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")[["3"]]
    probs <- probs[,"3"]
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    phe <- iron$pheno[,2,drop=FALSE]

    ind <- c(1:50, 101:150)
    expected <- scan1blup(probs[ind,], phe[ind,,drop=FALSE], kinship[ind,ind], addcovar=X[ind])
    expect_equal(scan1blup(probs[ind,], phe, kinship, addcovar=X), expected)
    expect_equal(scan1blup(probs, phe[ind,,drop=FALSE], kinship, addcovar=X), expected)
    expect_equal(scan1blup(probs, phe, kinship[ind,ind], addcovar=X), expected)
    expect_equal(scan1blup(probs, phe, kinship, addcovar=X[ind]), expected)

    # repeat with no kinship
    expected <- scan1blup(probs[ind,], phe[ind,,drop=FALSE], addcovar=X[ind])
    expect_equal(scan1blup(probs[ind,], phe, addcovar=X), expected)
    expect_equal(scan1blup(probs, phe[ind,,drop=FALSE], addcovar=X), expected)
    expect_equal(scan1blup(probs, phe, addcovar=X[ind]), expected)

})

context("BLUP estimates of QTL effects")

library(qtl2geno)
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
phe <- iron$pheno[,1,drop=FALSE]

test_that("scan1blup works with no kinship matrix", {

    pr <- calc_genoprob(iron[,"16"], step=0)

    blup <- scan1blup(pr, phe)
    blup_se <- scan1blup(pr, phe, se=TRUE)
    expect_equal(blup$coef, blup_se$coef)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, preserve_intercept=TRUE)
    expect_equal(blup$coef, blup_int$coef[,1:3] + blup_int$coef[,4])

    # brute force BLUPs
    for(i in 1:dim(pr$probs[[1]])[[3]]) {
        Z <- pr$probs[[1]][,,i]
        k <- Z %*% t(Z)
        ke <- decomp_kinship(k)
        keval <- ke$values
        kevec <- ke$vectors
        X <- cbind(rep(1, nrow(phe)))
        yp <- kevec %*% phe
        Xp <- kevec %*% X
        lmmfit <- Rcpp_fitLMM(keval, yp, Xp, tol=1e-12)

        beta <- lmmfit$beta
        expect_equal(blup_int$coef[i,4], beta, tolerance=1e-6)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (yp - Xp %*% beta) )
        names(u) <- colnames(blup$coef)
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-5)
    }

    # repeat with an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, addcovar=sex)
    blup_se <- scan1blup(pr, phe, addcovar=sex, se=TRUE)
    expect_equal(blup$coef, blup_se$coef)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, addcovar=sex, preserve_intercept=TRUE)
    expect_equal(blup$coef, blup_int$coef[,c(1:3,5)] + cbind(blup_int$coef[,c(4,4,4)], 0))

    # brute force BLUPs
    for(i in 1:dim(pr$probs[[1]])[[3]]) {
        Z <- pr$probs[[1]][,,i]
        k <- Z %*% t(Z)
        ke <- decomp_kinship(k)
        keval <- ke$values
        kevec <- ke$vectors
        X <- cbind(1, sex)
        yp <- kevec %*% phe
        Xp <- kevec %*% X
        lmmfit <- Rcpp_fitLMM(keval, yp, Xp, tol=1e-12)

        beta <- lmmfit$beta
        names(beta) <- c("intercept", "ac1")
        expect_equal(blup_int$coef[i,4:5], beta, tolerance=1e-6)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (yp - Xp %*% beta) )
        names(u) <- colnames(blup$coef)[1:3]
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-6)
    }

})

test_that("scan1blup works with kinship matrix", {

    pr <- calc_genoprob(iron, step=0)
    K <- calc_kinship(pr[,c(1:15,17:19,"X")])
    pr <- pr[,"16"]

    # sex as an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, K, sex)
    blup_se <- scan1blup(pr, phe, K, sex, se=TRUE)
    expect_equal(blup$coef, blup_se$coef)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, K, sex, preserve_intercept=TRUE)
    expect_equal(blup$coef, blup_int$coef[,c(1:3,5)] + cbind(blup_int$coef[,c(4,4,4)], 0))

    # brute force BLUPs
    Ke <- decomp_kinship(K)
    Keval <- Ke$values
    Kevec <- Ke$vectors
    yp <- Kevec %*% phe
    X <- cbind(1, sex)
    Xp <- Kevec %*% X
    lmmfit <- Rcpp_fitLMM(Keval, yp, Xp, tol=1e-12)
    hsq <- lmmfit$hsq
    wts <- 1/sqrt(hsq * Keval + (1-hsq))
    yp <- wts * yp
    Xp <- wts * Xp
    pradj <- pr$probs[[1]]
    for(i in 1:dim(pradj)[[3]])
        pradj[,,i] <- Kevec %*% pradj[,,i] * wts

    for(i in 1:dim(pradj)[[3]]) {
        Z <- pradj[,,i]
        k <- Z %*% t(Z)
        ke <- decomp_kinship(k)
        keval <- ke$values
        kevec <- ke$vectors
        ypp <- kevec %*% yp
        Xpp <- kevec %*% Xp
        lmmfit <- Rcpp_fitLMM(keval, ypp, Xpp, tol=1e-12)

        beta <- lmmfit$beta
        names(beta) <- c("intercept", "ac1")
        expect_equal(blup_int$coef[i,4:5], beta, tolerance=1e-3)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (ypp - Xpp %*% beta) )
        names(u) <- colnames(blup$coef)[1:3]
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-3)
    }

})

test_that("scan1blup works with kinship matrix on another chromosome", {

    pr <- calc_genoprob(iron, step=0)
    K <- calc_kinship(pr[,c(1:11,13:19,"X")])
    pr <- pr[,"12"]

    # sex as an additive covariate
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)

    blup <- scan1blup(pr, phe, K, sex)
    blup_se <- scan1blup(pr, phe, K, sex, se=TRUE)
    expect_equal(blup$coef, blup_se$coef)

    # cf preserve intercept vs not
    blup_int <- scan1blup(pr, phe, K, sex, preserve_intercept=TRUE)
    expect_equal(blup$coef, blup_int$coef[,c(1:3,5)] + cbind(blup_int$coef[,c(4,4,4)], 0))

    # brute force BLUPs
    Ke <- decomp_kinship(K)
    Keval <- Ke$values
    Kevec <- Ke$vectors
    yp <- Kevec %*% phe
    X <- cbind(1, sex)
    Xp <- Kevec %*% X
    lmmfit <- Rcpp_fitLMM(Keval, yp, Xp, tol=1e-12)
    hsq <- lmmfit$hsq
    wts <- 1/sqrt(hsq * Keval + (1-hsq))
    yp <- wts * yp
    Xp <- wts * Xp
    pradj <- pr$probs[[1]]
    for(i in 1:dim(pradj)[[3]])
        pradj[,,i] <- Kevec %*% pradj[,,i] * wts

    for(i in 1:dim(pradj)[[3]]) {
        Z <- pradj[,,i]
        k <- Z %*% t(Z)
        ke <- decomp_kinship(k)
        keval <- ke$values
        kevec <- ke$vectors
        ypp <- kevec %*% yp
        Xpp <- kevec %*% Xp
        lmmfit <- Rcpp_fitLMM(keval, ypp, Xpp, tol=1e-12)

        beta <- lmmfit$beta
        names(beta) <- c("intercept", "ac1")
        expect_equal(blup_int$coef[i,4:5], beta, tolerance=1e-3)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (ypp - Xpp %*% beta) )
        names(u) <- colnames(blup$coef)[1:3]
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-3)
    }

})

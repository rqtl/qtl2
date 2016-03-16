context("BLUP estimates of QTL effects")

test_that("scan1blup works with no kinship matrix", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    pr <- calc_genoprob(iron[,"16"], step=0)
    phe <- iron$pheno[,1,drop=FALSE]

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
        expect_equal(blup_int$coef[i,4], beta)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (yp - Xp %*% beta) )
        names(u) <- colnames(blup$coef)
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-6)
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
        expect_equal(blup_int$coef[i,4:5], beta)

        hsq <- lmmfit$hsq
        wts <- hsq * keval + (1-hsq)
        u <- as.numeric( t(kevec %*% Z) %*% diag(hsq/wts) %*%  (yp - Xp %*% beta) )
        names(u) <- colnames(blup$coef)[1:3]
        expect_equal(blup_int$coef[i, 1:3], u, tolerance=1e-6)
    }

})

context("LMM genome scan by scan1_lmm")


test_that("scan1_lmm with intercross, vs ported lmmlite code", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)

    out_reml <- scan1_lmm(probs, iron$pheno, kinship, reml=TRUE)
    out_ml <- scan1_lmm(probs, iron$pheno, kinship, reml=FALSE)

    # "by hand" calculation
    y <- iron$pheno
    X <- cbind(rep(1, nrow(iron$pheno)))
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% X

    byhand1_reml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=TRUE, tol=1e-12)
    byhand2_reml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=TRUE, tol=1e-12)
    byhand1_ml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=FALSE, tol=1e-12)
    byhand2_ml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=FALSE, tol=1e-12)

    # hsq the same?
    expect_equal(as.numeric(attr(out_reml, "hsq")[1,]),
                 c(byhand1_reml$hsq, byhand2_reml$hsq))
    expect_equal(as.numeric(attr(out_ml, "hsq")[1,]),
                 c(byhand1_ml$hsq, byhand2_ml$hsq))

    # compare chromosome 1 LOD scores
    d <- dim(probs[[1]])[3]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, d)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(X, probs[[1]][,-1,i])
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2_ml$loglik)/log(10)

    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[1:d,1], lod_reml1)
    expect_equal(out_reml[1:d,2], lod_reml2)
    expect_equal(out_ml[1:d,1], lod_ml1)
    expect_equal(out_ml[1:d,2], lod_ml2)

})

test_that("scan1_lmm with intercross with X covariates for null", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)

    out_reml <- scan1_lmm(probs, iron$pheno, kinship, Xcovar=Xc, reml=TRUE)
    out_ml <- scan1_lmm(probs, iron$pheno, kinship, Xcovar=Xc, reml=FALSE)

    # "by hand" calculation
    y <- iron$pheno
    X <- cbind(rep(1, nrow(iron$pheno)))
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% X
    Xcp <- Ke$vectors %*% Xc

    byhand1_reml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp), reml=TRUE, tol=1e-12)
    byhand2_reml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp), reml=TRUE, tol=1e-12)
    byhand1_ml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp), reml=FALSE, tol=1e-12)
    byhand2_ml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp), reml=FALSE, tol=1e-12)

    # hsq the same?
    expect_equal(as.numeric(attr(out_reml, "hsq")[2,]),
                 c(byhand1_reml$hsq, byhand2_reml$hsq), tolerance=1e-6)
    expect_equal(as.numeric(attr(out_ml, "hsq")[2,]),
                 c(byhand1_ml$hsq, byhand2_ml$hsq))

    # compare chromosome X LOD scores
    d <- dim(probs[["X"]])[3]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, d)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(1, probs[["X"]][,-1,i])
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2_ml$loglik)/log(10)

    index <- nrow(out_reml) - rev(1:d) + 1
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[index,1], lod_reml1)
    expect_equal(out_reml[index,2], lod_reml2, tolerance=1e-6)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

})


test_that("scan1_lmm with intercross with an additive covariate", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_ml <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc, reml=FALSE, tol=1e-12)

    # "by hand" calculation
    y <- iron$pheno
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% cbind(1, X)
    Xcp <- Ke$vectors %*% Xc

    # autosome null
    byhand1A_reml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=TRUE, tol=1e-12)
    byhand2A_reml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=TRUE, tol=1e-12)
    byhand1A_ml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=FALSE, tol=1e-12)
    byhand2A_ml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=FALSE, tol=1e-12)

    expect_equal(as.numeric(attr(out_reml, "hsq")[1,]),
                 c(byhand1A_reml$hsq, byhand2A_reml$hsq))
    expect_equal(as.numeric(attr(out_ml, "hsq")[1,]),
                 c(byhand1A_ml$hsq, byhand2A_ml$hsq))

    # X chr null
    byhand1X_reml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp[,-1]), reml=TRUE, tol=1e-12)
    byhand2X_reml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp[,-1]), reml=TRUE, tol=1e-12)
    byhand1X_ml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp[,-1]), reml=FALSE, tol=1e-12)
    byhand2X_ml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp[,-1]), reml=FALSE, tol=1e-12)

    # hsq the same?
    expect_equal(as.numeric(attr(out_reml, "hsq")[2,]),
                 c(byhand1X_reml$hsq, byhand2X_reml$hsq), tolerance=1e-6)
    expect_equal(as.numeric(attr(out_ml, "hsq")[2,]),
                 c(byhand1X_ml$hsq, byhand2X_ml$hsq))


    # compare chromosome 2 LOD scores
    d <- dim(probs[["2"]])[3]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, d)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(1, X, probs[["2"]][,-1,i])
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1A_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2A_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1A_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2A_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1A_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2A_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1A_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2A_ml$loglik)/log(10)

    index <- dim(probs[["1"]])[3] + 1:dim(probs[["2"]])[3]
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[index,1], lod_reml1)
    expect_equal(out_reml[index,2], lod_reml2)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

    # compare chromosome X LOD scores
    d <- dim(probs[["X"]])[3]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, d)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(1, probs[["X"]][,-1,i])
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1X_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2X_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1X_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2X_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1X_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2X_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1X_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2X_ml$loglik)/log(10)

    index <- nrow(out_reml) - rev(1:d) + 1
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    ## FIX_ME
    ## REML not yet working on X chromosome, when (X, probs) is not full rank
#    expect_equal(out_reml[index,1], lod_reml1)
#    expect_equal(out_reml[index,2], lod_reml2)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

})


test_that("scan1_lmm with intercross with an interactive covariate", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=2.5, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                          Xcovar=Xc, reml=TRUE, tol=1e-12, intcovar_method="lowmem")
    out_ml <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                        Xcovar=Xc, reml=FALSE, tol=1e-12, intcovar_method="lowmem")
    out_reml_himem <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                          Xcovar=Xc, reml=TRUE, tol=1e-12, intcovar_method="highmem")
    out_ml_himem <- scan1_lmm(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                        Xcovar=Xc, reml=FALSE, tol=1e-12, intcovar_method="highmem")

    # same result using "highmem" and "lowmem" methods
    expect_equal(out_reml_himem, out_reml)
    expect_equal(out_ml_himem, out_ml)

    # "by hand" calculation
    y <- iron$pheno
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% cbind(1, X)
    Xcp <- Ke$vectors %*% Xc

    # autosome null (same as w/o interactive covariate)
    byhand1A_reml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=TRUE, tol=1e-12)
    byhand2A_reml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=TRUE, tol=1e-12)
    byhand1A_ml <- Rcpp_fitLMM(Ke$values, yp[,1], Xp, reml=FALSE, tol=1e-12)
    byhand2A_ml <- Rcpp_fitLMM(Ke$values, yp[,2], Xp, reml=FALSE, tol=1e-12)

    expect_equal(as.numeric(attr(out_reml, "hsq")[1,]),
                 c(byhand1A_reml$hsq, byhand2A_reml$hsq))
    expect_equal(as.numeric(attr(out_ml, "hsq")[1,]),
                 c(byhand1A_ml$hsq, byhand2A_ml$hsq))

    # X chr null (same as w/o interactive covariate)
    byhand1X_reml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp[,-1]), reml=TRUE, tol=1e-12)
    byhand2X_reml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp[,-1]), reml=TRUE, tol=1e-12)
    byhand1X_ml <- Rcpp_fitLMM(Ke$values, yp[,1], cbind(Xp, Xcp[,-1]), reml=FALSE, tol=1e-12)
    byhand2X_ml <- Rcpp_fitLMM(Ke$values, yp[,2], cbind(Xp, Xcp[,-1]), reml=FALSE, tol=1e-12)

    # hsq the same?
    expect_equal(as.numeric(attr(out_reml, "hsq")[2,]),
                 c(byhand1X_reml$hsq, byhand2X_reml$hsq), tolerance=1e-6)
    expect_equal(as.numeric(attr(out_ml, "hsq")[2,]),
                 c(byhand1X_ml$hsq, byhand2X_ml$hsq))


    # compare chromosome 4 LOD scores
    npos <- sapply(probs, function(a) dim(a)[[3]])
    d <- npos["4"]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, d)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(1, X, probs[["4"]][,-1,i], probs[["4"]][,-1,i]*X)
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1A_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2A_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1A_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2A_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1A_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2A_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1A_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2A_ml$loglik)/log(10)

    index <- sum(npos[1:3]) + 1:npos[4]
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[index,1], lod_reml1)
    expect_equal(out_reml[index,2], lod_reml2)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

    # compare chromosome X LOD scores
    d <- dim(probs[["X"]])[3]
    loglik_reml1 <- loglik_reml2 <-
        loglik_ml1 <- loglik_ml2 <- rep(NA, 3)
    for(i in 1:d) {
        Xp <- Ke$vectors %*% cbind(1, X, probs[["X"]][,-1,i], probs[["X"]][,-1,i]*X)
        # calculate likelihoods using plain ML (not the residual log likelihood)
        loglik_reml1[i] <- Rcpp_calcLL(byhand1X_reml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_reml2[i] <- Rcpp_calcLL(byhand2X_reml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
        loglik_ml1[i] <- Rcpp_calcLL(byhand1X_ml$hsq, Ke$values, yp[,1], Xp, reml=FALSE)
        loglik_ml2[i] <- Rcpp_calcLL(byhand2X_ml$hsq, Ke$values, yp[,2], Xp, reml=FALSE)
    }
    lod_reml1 <- (loglik_reml1 - byhand1X_reml$loglik)/log(10)
    lod_reml2 <- (loglik_reml2 - byhand2X_reml$loglik)/log(10)
    lod_ml1 <- (loglik_ml1 - byhand1X_ml$loglik)/log(10)
    lod_ml2 <- (loglik_ml2 - byhand2X_ml$loglik)/log(10)

    index <- nrow(out_reml) - rev(1:d) + 1
    ## FIX ME
    ## FIX_ME
    ## Not yet working on X chromosome, when (X, probs) is not full rank
#    expect_equal(out_reml[index,1], lod_reml1)
#    expect_equal(out_reml[index,2], lod_reml2)
#    expect_equal(out_ml[index,1], lod_ml1)
#    expect_equal(out_ml[index,2], lod_ml2)

})

# look at LOCO

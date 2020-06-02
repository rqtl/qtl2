context("LMM genome scan by scan1 with kinship matrix")


test_that("scan1 with kinship with intercross, vs ported lmmlite code", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)

    out_reml <- scan1(probs, iron$pheno, kinship, reml=TRUE)
    out_ml <- scan1(probs, iron$pheno, kinship, reml=FALSE)

    # "by hand" calculation
    y <- iron$pheno
    X <- cbind(rep(1, nrow(iron$pheno)))
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% X
    # double the eigenvalues (== kinship matrix * 2)
    Ke$values <- Ke$values*2

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

    out_reml <- unclass(out_reml)
    out_ml <- unclass(out_ml)
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[1:d,1], lod_reml1)
    expect_equal(out_reml[1:d,2], lod_reml2)
    expect_equal(out_ml[1:d,1], lod_ml1)
    expect_equal(out_ml[1:d,2], lod_ml2)

})

test_that("scan1 with intercross with X covariates for null", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)

    out_reml <- scan1(probs, iron$pheno, kinship, Xcovar=Xc, reml=TRUE)
    out_ml <- scan1(probs, iron$pheno, kinship, Xcovar=Xc, reml=FALSE)

    # "by hand" calculation
    y <- iron$pheno
    X <- cbind(rep(1, nrow(iron$pheno)))
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% X
    Xcp <- Ke$vectors %*% Xc
    # double the eigenvalues (== kinship matrix * 2)
    Ke$values <- Ke$values*2

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
    out_reml <- unclass(out_reml)
    out_ml <- unclass(out_ml)
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    expect_equal(out_reml[index,1], lod_reml1, tolerance=1e-6)
    expect_equal(out_reml[index,2], lod_reml2, tolerance=1e-6)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

})


test_that("scan1 with kinship with intercross with an additive covariate", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc, reml=FALSE, tol=1e-12)

    # "by hand" calculation
    y <- iron$pheno
    Ke <- decomp_kinship(kinship) # eigen decomp
    yp <- Ke$vectors %*% y
    Xp <- Ke$vectors %*% cbind(1, X)
    Xcp <- Ke$vectors %*% Xc
    # double the eigenvalues (== kinship matrix * 2)
    Ke$values <- Ke$values*2

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
    out_reml <- unclass(out_reml)
    out_ml <- unclass(out_ml)
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
    out_reml <- unclass(out_reml)
    out_ml <- unclass(out_ml)
    dimnames(out_reml) <- dimnames(out_ml) <- NULL
    ## FIX_ME
    ## REML not yet working on X chromosome, when (X, probs) is not full rank
#    expect_equal(out_reml[index,1], lod_reml1)
#    expect_equal(out_reml[index,2], lod_reml2)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

})


test_that("scan1 with kinship with intercross with an interactive covariate", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12, intcovar_method="lowmem")
    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12, intcovar_method="lowmem")
    out_reml_himem <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                            Xcovar=Xc, reml=TRUE, tol=1e-12, intcovar_method="highmem")
    out_ml_himem <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
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
    # double the eigenvalues (== kinship matrix * 2)
    Ke$values <- Ke$values*2

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
    out_reml <- unclass(out_reml)
    out_ml <- unclass(out_ml)
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
    ## Not yet working on X chromosome, when (X, probs) is not full rank
#    expect_equal(out_reml[index,1], lod_reml1)
#    expect_equal(out_reml[index,2], lod_reml2)
#    expect_equal(out_ml[index,1], lod_ml1)
#    expect_equal(out_ml[index,2], lod_ml2)

})

test_that("scan1 with kinship works with LOCO, additive covariates", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12)

    y <- iron$pheno
    Ke <- decomp_kinship(kinship) # eigen decomp
    # double the eigenvalues (== kinship matrix * 2)
    Ke <- lapply(Ke, function(a) { a$values <- 2*a$values; a})

    # compare chromosomes 1, 6, 9, 18
    chrs <- paste(c(1,6,9,18))
    npos <- sapply(probs, function(a) dim(a)[[3]])

    for(chr in chrs) {
        nchr <- which(names(npos) == chr)
        d <- npos[chr]

        yp <- Ke[[chr]]$vectors %*% y
        Xp <- Ke[[chr]]$vectors %*% cbind(1, X)

        # autosome null
        byhand1_reml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,1], Xp, reml=TRUE, tol=1e-12)
        byhand2_reml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,2], Xp, reml=TRUE, tol=1e-12)
        byhand1_ml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,1], Xp, reml=FALSE, tol=1e-12)
        byhand2_ml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,2], Xp, reml=FALSE, tol=1e-12)

        expect_equal(as.numeric(attr(out_reml, "hsq")[nchr,]),
                     c(byhand1_reml$hsq, byhand2_reml$hsq), tolerance=1e-5)
        expect_equal(as.numeric(attr(out_ml, "hsq")[nchr,]),
                     c(byhand1_ml$hsq, byhand2_ml$hsq), tolerance=1e-6)

        # chromosome scan
        loglik_reml1 <- loglik_reml2 <-
            loglik_ml1 <- loglik_ml2 <- rep(NA, d)
        for(i in 1:d) {
            Xp <- Ke[[chr]]$vectors %*% cbind(1, X, probs[[chr]][,-1,i])
            # calculate likelihoods using plain ML (not the residual log likelihood)
            loglik_reml1[i] <- Rcpp_calcLL(byhand1_reml$hsq, Ke[[chr]]$values, yp[,1], Xp, reml=FALSE)
            loglik_reml2[i] <- Rcpp_calcLL(byhand2_reml$hsq, Ke[[chr]]$values, yp[,2], Xp, reml=FALSE)
            loglik_ml1[i] <- Rcpp_calcLL(byhand1_ml$hsq, Ke[[chr]]$values, yp[,1], Xp, reml=FALSE)
            loglik_ml2[i] <- Rcpp_calcLL(byhand2_ml$hsq, Ke[[chr]]$values, yp[,2], Xp, reml=FALSE)
        }
        lod_reml1 <- (loglik_reml1 - byhand1_reml$loglik)/log(10)
        lod_reml2 <- (loglik_reml2 - byhand2_reml$loglik)/log(10)
        lod_ml1 <- (loglik_ml1 - byhand1_ml$loglik)/log(10)
        lod_ml2 <- (loglik_ml2 - byhand2_ml$loglik)/log(10)

        if(nchr > 1) index <- sum(npos[1:(nchr-1)]) + 1:d
        else index <- 1:d
        out_reml <- unclass(out_reml)
        out_ml <- unclass(out_ml)
        dimnames(out_reml) <- dimnames(out_ml) <- NULL
        expect_equal(out_reml[index,1], lod_reml1)
        expect_equal(out_reml[index,2], lod_reml2, tolerance=1e-5)
        expect_equal(out_ml[index,1], lod_ml1)
        expect_equal(out_ml[index,2], lod_ml2)
    }

})

test_that("scan1 with kinship works with LOCO, interactive covariates", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12)


    y <- iron$pheno
    Ke <- decomp_kinship(kinship) # eigen decomp
    # double the eigenvalues (== kinship matrix * 2)
    Ke <- lapply(Ke, function(a) { a$values <- 2*a$values; a})

    # compare chromosomes 1, 6, 9, 18
    chrs <- paste(c(1,6,9,18))
    npos <- sapply(probs, function(a) dim(a)[[3]])

    for(chr in chrs) {
        nchr <- which(names(npos) == chr)
        d <- npos[chr]

        yp <- Ke[[chr]]$vectors %*% y
        Xp <- Ke[[chr]]$vectors %*% cbind(1, X)

        # autosome null (same as w/o interactive covariate)
        byhand1_reml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,1], Xp, reml=TRUE, tol=1e-12)
        byhand2_reml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,2], Xp, reml=TRUE, tol=1e-12)
        byhand1_ml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,1], Xp, reml=FALSE, tol=1e-12)
        byhand2_ml <- Rcpp_fitLMM(Ke[[chr]]$values, yp[,2], Xp, reml=FALSE, tol=1e-12)

        expect_equal(as.numeric(attr(out_reml, "hsq")[nchr,]),
                     c(byhand1_reml$hsq, byhand2_reml$hsq), tolerance=1e-5)
        expect_equal(as.numeric(attr(out_ml, "hsq")[nchr,]),
                     c(byhand1_ml$hsq, byhand2_ml$hsq), tolerance=1e-6)

        # chromosome scan
        loglik_reml1 <- loglik_reml2 <-
            loglik_ml1 <- loglik_ml2 <- rep(NA, d)
        for(i in 1:d) {
            Xp <- Ke[[chr]]$vectors %*% cbind(1, X, probs[[chr]][,-1,i], probs[[chr]][,-1,i]*X)
            # calculate likelihoods using plain ML (not the residual log likelihood)
            loglik_reml1[i] <- Rcpp_calcLL(byhand1_reml$hsq, Ke[[chr]]$values, yp[,1], Xp, reml=FALSE)
            loglik_reml2[i] <- Rcpp_calcLL(byhand2_reml$hsq, Ke[[chr]]$values, yp[,2], Xp, reml=FALSE)
            loglik_ml1[i] <- Rcpp_calcLL(byhand1_ml$hsq, Ke[[chr]]$values, yp[,1], Xp, reml=FALSE)
            loglik_ml2[i] <- Rcpp_calcLL(byhand2_ml$hsq, Ke[[chr]]$values, yp[,2], Xp, reml=FALSE)
        }
        lod_reml1 <- (loglik_reml1 - byhand1_reml$loglik)/log(10)
        lod_reml2 <- (loglik_reml2 - byhand2_reml$loglik)/log(10)
        lod_ml1 <- (loglik_ml1 - byhand1_ml$loglik)/log(10)
        lod_ml2 <- (loglik_ml2 - byhand2_ml$loglik)/log(10)

        if(nchr > 1) index <- sum(npos[1:(nchr-1)]) + 1:d
        else index <- 1:d
        out_reml <- unclass(out_reml)
        out_ml <- unclass(out_ml)
        dimnames(out_reml) <- dimnames(out_ml) <- NULL
        expect_equal(out_reml[index,1], lod_reml1)
        expect_equal(out_reml[index,2], lod_reml2, tolerance=1e-6)
        expect_equal(out_ml[index,1], lod_ml1)
        expect_equal(out_ml[index,2], lod_ml2)
    }


})


test_that("scan1 with kinship works with multicore", {

    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_reml_4core <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                            Xcovar=Xc, reml=TRUE, tol=1e-12, cores=2)
    expect_equal(out_reml, out_reml_4core)


    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12)
    out_ml_4core <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                          Xcovar=Xc, reml=FALSE, tol=1e-12, cores=2)
    expect_equal(out_ml, out_ml_4core)

})


test_that("scan1 with kinship LOD results invariant to change in scale to pheno and covar", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12)
    out_reml_scale <- scan1(probs, iron$pheno/100, kinship, addcovar=X*2, intcovar=X*2,
                            Xcovar=Xc*2, reml=TRUE, tol=1e-12)
    expect_equal(out_reml, out_reml_scale, tol=1e-6)


    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12)
    out_ml_scale <- scan1(probs, iron$pheno/100, kinship, addcovar=X*4, intcovar=X*4,
                          Xcovar=Xc*4, reml=FALSE, tol=1e-12)
    expect_equal(out_ml, out_ml_scale, tol=1e-6)

})

test_that("scan1 deals with mismatching individuals", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    ind <- c(1:50, 101:150)
    subK <- lapply(kinship, "[", ind, ind)
    expected <- scan1(probs[ind,], iron$pheno[ind,,drop=FALSE], subK, addcovar=X[ind], intcovar=X[ind],
                      Xcovar=Xc[ind,], reml=TRUE, tol=1e-12)
    expect_equal(scan1(probs[ind,], iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12), expected)
    expect_equal(scan1(probs, iron$pheno[ind,], kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12), expected)
    expect_equal(scan1(probs, iron$pheno, subK, addcovar=X, intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12), expected)
    expect_equal(scan1(probs, iron$pheno, kinship, addcovar=X[ind], intcovar=X,
                      Xcovar=Xc, reml=TRUE, tol=1e-12), expected)
    expect_equal(scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X[ind],
                      Xcovar=Xc, reml=TRUE, tol=1e-12), expected)
    expect_equal(scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                      Xcovar=Xc[ind,], reml=TRUE, tol=1e-12), expected)

})


test_that("scan1 with weights and kinship", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    chr <- c("7", "12", "X")
    probs <- probs[,chr]
    kinship <- kinship[chr]

    RNGkind("Mersenne-Twister") # make sure we're using the standard RNG
    set.seed(28915967)
    weights <- stats::setNames(sample(1:10, n_ind(iron), replace=TRUE), ind_ids(iron))

    out_reml <- scan1(probs, iron$pheno, kinship, addcovar=X,
                      Xcovar=Xc, weights=weights, reml=TRUE, tol=1e-12)

    ### results via regress library
    # library(regress)
    # hsq <- matrix(nrow=3, ncol=2)
    # dimnames(hsq) <- list(c("7", "12", "X"), c("liver", "spleen"))
    # for(i in 1:2) {
    #     for(j in 1:3) {
    #         k <- kinship[[j]]*2
    #         if(j==3) co <- Xc else co <- X
    #         out_regress <- regress(iron$pheno[,i] ~ co, ~ k + diag(1/weights), identity=FALSE, tol=1e-8)
    #         sig <- out_regress$sigma
    #         hsq[j,i] <- sig[1]/sum(sig)
    #    }
    # }
    hsq <- structure(c(0.281878438466624, 0.353558919335662, 0.354315217715963,
                       0.20945078570429, 0.217671918799751, 0.201682575651958),
                     .Dim = 3:2, .Dimnames = list(c("7", "12", "X"), c("liver", "spleen")))

    expect_equal(attr(out_reml, "hsq"), hsq, tol=1e-5)

    ### calculate lod scors from scan
    # hsq <- attr(out_reml, "hsq")
    # result <- vector("list", 6)
    # for(i in 1:2) {
    #     for(j in 1:3) {
    #         # cholesky decomp of kinhip matrix
    #         v <- hsq[j,i]*2*kinship[[j]] + (1-hsq[j,i])*diag(1/weights)
    #         d <- solve(t(chol( v )))
    #         # pre-multiply phenotype and covar by t(d)
    #         y <- d %*% iron$pheno[,i]
    #         if(j==3) co <- Xc else co <- X
    #         co <- d %*% cbind(1, co)
    #
    #         # log lik under null
    #         r <- lm(y ~ -1 + co)$resid
    #         llik0 <- sum(dnorm(r, 0, sqrt(mean(r^2)), log=TRUE))
    #
    #         llik1 <- apply(probs[[j]], 3, function(z) {
    #             thisX <- cbind(co, d %*% z[,-1,drop=FALSE])
    #             r <- lm(y ~ -1 + thisX)$resid
    #             sum(dnorm(r, 0, sqrt(mean(r^2)), log=TRUE)) })
    #
    #         result[[(i-1)*3+j]] <- (llik1 - llik0)/log(10)
    #
    #     }
    # }
    # result <- rbind(cbind(result[[1]], result[[4]]),
    #                 cbind(result[[2]], result[[5]]),
    #                 cbind(result[[3]], result[[6]]))
    # dimnames(result) <- dimnames(out_reml)
    # class(result) <- class(out_reml)

    result <- structure(c(0.765379546370331, 1.33697830149196, 2.05141217249269,
                          2.82748389184348, 3.57825928288821, 4.11482576884623,
                          4.18208695088672, 4.5010472223741, 4.77640847444171,
                          4.98121670508953, 5.07224467963965, 5.28435499128555,
                          5.10772158608095, 3.83479814287182, 3.89931966500793,
                          4.24832760135268, 4.18515409936768, 4.61791967696627,
                          4.70308738505572, 4.58758777783123, 4.88425659502959,
                          5.39302146958257, 5.86125538183219, 6.27973204425698,
                          6.62245498126028, 6.84349072896544, 6.91418168584823,
                          0.492648140597935, 0.483021280268028, 0.465046322957041,
                          0.435712712087884, 0.391979816567206, 0.332013881699469,
                          0.257330299826467, 0.174985513110401, 0.0975448556296765,
                          0.0390292238934769, 0.0087954970835623, 0.0081071126160765,
                          0.0315063370227757, 0.0706853835370534, 0.117761588763882,
                          0.166806834360401, 0.180329873353274, 0.146302909441694,
                          0.161137797908763, 0.192675853649471, 0.243358929100302,
                          0.313640413846475, 0.401241925760445, 0.501210880937258,
                          0.606957139785506, 0.711817493234455, 0.810381904299065,
                          0.899112752704855, 0.976303366208711, 1.00118015673376,
                          0.68879116519848, 1.13593060562675, 1.59214070065802,
                          1.93410070637396, 2.08715532129677, 2.08161469298463,
                          2.11418090796465, 2.24129850312586, 2.28750006399176,
                          2.24152015532401, 2.15435105671462, 2.16059659637594,
                          1.40659065214545, 0.389546568631669, 0.392567989473253,
                          0.432129729053773, 0.445384879284572, 0.468919719042007,
                          0.561846602500791, 0.61487493505232, 0.680306211651546,
                          0.828892839652596, 1.03747954481922, 1.3178288863555,
                          1.64752359905175, 1.97570834069646, 2.25872692788406,
                          1.76552836513385, 1.9624748694232, 2.179534029569,
                          2.40709075366897, 2.62706355895772, 2.81242765604838,
                          2.93166656075001, 2.95963946310072, 2.89030386259124,
                          2.74110806468527, 2.54394356128059, 2.33040574729971,
                          2.12258054254516, 1.93206224497564, 1.76308130810039,
                          1.61583945469733, 1.57830829020771, 1.91298146944818,
                          1.89937590028851, 1.86320092206996, 1.79734091344691,
                          1.69690361501555, 1.56240479342725, 1.4020247631517,
                          1.23069582035597, 1.06561266027339, 0.920732232061234,
                          0.803561188099648, 0.715279174632663, 0.690066084958058),
                        .Dim = c(57L, 2L),
                        .Dimnames = list(c("D7Mit74", "c7.loc4", "c7.loc6", "c7.loc9",
                                           "c7.loc11", "D7Mit25", "c7.loc14", "c7.loc16",
                                           "c7.loc19", "c7.loc21", "D7Nds5", "c7.loc24",
                                           "c7.loc26", "D7mit30", "c7.loc29", "c7.loc31",
                                           "D7Mit31", "c7.loc34", "c7.loc36", "D7Mit17",
                                           "c7.loc39", "c7.loc41", "c7.loc44", "c7.loc46",
                                           "c7.loc49", "c7.loc51", "D7Mit71", "D12Mit88",
                                           "c12.loc22", "c12.loc25", "c12.loc27", "c12.loc30",
                                           "c12.loc32", "c12.loc35", "c12.loc37", "c12.loc40",
                                           "c12.loc42", "c12.loc45", "c12.loc47", "c12.loc50",
                                           "c12.loc52", "c12.loc55", "c12.loc57", "D12Mit134",
                                           "DXMit16", "cX.loc32", "cX.loc34", "cX.loc37",
                                           "cX.loc40", "cX.loc42", "cX.loc44", "cX.loc47",
                                           "cX.loc50", "cX.loc52", "cX.loc54",
                                           "cX.loc57", "DXMit186"), c("liver", "spleen")),
                        class = c("scan1", "matrix"))

    for(at in c("hsq", "sample_size"))
        attr(result, at) <- attr(out_reml, at)
    expect_equal(out_reml, result, tol=1e-6)

})

test_that("scan1 works with hsq specified", {

    if(isnt_karl()) skip("This test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    kinship_loco <- calc_kinship(probs, "loco")

    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    # no covariates, single kinship
    out <- scan1(probs, iron$pheno, kinship)
    expect_equal(scan1(probs, iron$pheno, kinship, hsq=attr(out, "hsq")), out)

    # no covariates, loco kinship
    out <- scan1(probs, iron$pheno, kinship_loco)
    expect_equal(scan1(probs, iron$pheno, kinship_loco, hsq=attr(out, "hsq")), out)


    # X covariates, single kinship
    out <- scan1(probs, iron$pheno, kinship, Xcovar=Xc)
    expect_equal(scan1(probs, iron$pheno, kinship, Xcovar=Xc, hsq=attr(out, "hsq")), out)

    # X covariates, loco kinship
    out <- scan1(probs, iron$pheno, kinship_loco, Xcovar=Xc)
    expect_equal(scan1(probs, iron$pheno, kinship_loco, Xcovar=Xc, hsq=attr(out, "hsq")), out)


    # covariates, single kinship
    out <- scan1(probs, iron$pheno, kinship, addcovar=X)
    expect_equal(scan1(probs, iron$pheno, kinship, addcovar=X, hsq=attr(out, "hsq")), out)

    # covariates, loco kinship
    out <- scan1(probs, iron$pheno, kinship_loco, addcovar=X)
    expect_equal(scan1(probs, iron$pheno, kinship_loco, addcovar=X, hsq=attr(out, "hsq")), out)



    # covariates + X covar, single kinship
    out <- scan1(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc)
    expect_equal(scan1(probs, iron$pheno, kinship, addcovar=X, Xcovar=Xc, hsq=attr(out, "hsq")), out)

    # covariates + X covar, loco kinship
    out <- scan1(probs, iron$pheno, kinship_loco, addcovar=X, Xcovar=Xc)
    expect_equal(scan1(probs, iron$pheno, kinship_loco, addcovar=X, Xcovar=Xc, hsq=attr(out, "hsq")), out)

})

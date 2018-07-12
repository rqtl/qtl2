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
    expect_equal(out_reml[index,1], lod_reml1)
    expect_equal(out_reml[index,2], lod_reml2, tolerance=1e-6)
    expect_equal(out_ml[index,1], lod_ml1)
    expect_equal(out_ml[index,2], lod_ml2)

})


test_that("scan1 with kinship with intercross with an additive covariate", {

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
                            Xcovar=Xc, reml=TRUE, tol=1e-12, cores=4)
    expect_equal(out_reml, out_reml_4core)


    out_ml <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                    Xcovar=Xc, reml=FALSE, tol=1e-12)
    out_ml_4core <- scan1(probs, iron$pheno, kinship, addcovar=X, intcovar=X,
                          Xcovar=Xc, reml=FALSE, tol=1e-12, cores=4)
    expect_equal(out_ml, out_ml_4core)

})


test_that("scan1 with kinship LOD results invariant to change in scale to pheno and covar", {
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
    hsq <- structure(c(0.294000449087691, 0.386627569931694, 0.378245886332897,
                            0.175506814400581, 0.168504256116417, 0.162320444038263),
                          .Dim = c(3L, 2L), .Dimnames = list(c("7", "12", "X"), c("liver", "spleen")))
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

    result <- structure(c(0.76665906010563, 1.24382158309899, 1.81224411849879,
                          2.40917983502838, 2.9681402867112, 3.34660404416498, 3.40242177199974,
                          3.64993252963335, 3.82831170740312, 3.92278298641073, 3.93646354712277,
                          4.22771457688576, 4.31091438174264, 2.87038565955201, 2.92332212813056,
                          3.38043927872609, 3.40984534550206, 3.69585404453244, 3.83202186149427,
                          3.81751975419027, 4.13854424550092, 4.75476343452317, 5.40346726961252,
                          6.01943949217954, 6.50227740212975, 6.77330850173292, 6.82326722876252,
                          0.0337635407404249, 0.0295073664235682, 0.0248963840282653, 0.0201302579703893,
                          0.0155434921715859, 0.0115999393011996, 0.00880711811596259,
                          0.00753065584579277, 0.00778249200204616, 0.00914745375294169,
                          0.0109485496086686, 0.01254381194159, 0.0135420319954184, 0.0138371412928176,
                          0.0135187847307913, 0.0127601407700048, 0.0124937220029721, 0.34084397693496,
                          0.229366404351028, 0.140459622318806, 0.092548361857524, 0.104314974037626,
                          0.188035109538113, 0.342842577101976, 0.552805457762931, 0.7921888405597,
                          1.03460774402442, 1.26024748384826, 1.45822405357864, 1.52189801223226,
                          0.910960660359077, 1.4173925059434, 1.92335855116895, 2.24791673338013,
                          2.26522372565359, 2.07243617273106, 2.08641704560653, 2.1193512771495,
                          2.06631049680757, 1.9156254200745, 1.7498499850465, 1.5934400698397,
                          0.635554443348193, 0.290528704527623, 0.278241562971952, 0.136474241619778,
                          0.111524098169851, 0.0334455140979009, 0.00185536203197257, 0.017075001515959,
                          0.03595993711545, 0.106071668059052, 0.237852958970144, 0.437850802828232,
                          0.691690742231252, 0.970480263053637, 1.24748637774199, 1.2287845069906,
                          1.32089960813426, 1.41704748977607, 1.51223775484478, 1.59880158908385,
                          1.66680389654129, 1.70601473711322, 1.70938856968217, 1.67642327126315,
                          1.61395984486654, 1.53354837522464, 1.44715601460968, 1.36385997346139,
                          1.28876904549778, 1.22367673838179, 1.16832881654631, 1.15442476120788,
                          2.16884478394936, 2.19960890665743, 2.23654790432899, 2.27849702371477,
                          2.32300613699966, 2.36662725497101, 2.40582589433677, 2.43811505757396,
                          2.46266497024927, 2.48004435865592, 2.49149400716817, 2.49832496179339,
                          2.499876406515), .Dim = c(57L, 2L),
                        .Dimnames = list(c("D7Mit74",
                                           "c7.loc4", "c7.loc6", "c7.loc9", "c7.loc11", "D7Mit25", "c7.loc14",
                                           "c7.loc16", "c7.loc19", "c7.loc21", "D7Nds5", "c7.loc24", "c7.loc26",
                                           "D7mit30", "c7.loc29", "c7.loc31", "D7Mit31", "c7.loc34", "c7.loc36",
                                           "D7Mit17", "c7.loc39", "c7.loc41", "c7.loc44", "c7.loc46", "c7.loc49",
                                           "c7.loc51", "D7Mit71", "D12Mit88", "c12.loc22", "c12.loc25",
                                           "c12.loc27", "c12.loc30", "c12.loc32", "c12.loc35", "c12.loc37",
                                           "c12.loc40", "c12.loc42", "c12.loc45", "c12.loc47", "c12.loc50",
                                           "c12.loc52", "c12.loc55", "c12.loc57", "D12Mit134", "DXMit16",
                                           "cX.loc32", "cX.loc34", "cX.loc37", "cX.loc40", "cX.loc42", "cX.loc44",
                                           "cX.loc47", "cX.loc50", "cX.loc52", "cX.loc54", "cX.loc57", "DXMit186"),
                                         c("liver", "spleen")), class = c("scan1", "matrix"))

    for(at in c("hsq", "sample_size"))
        attr(result, at) <- attr(out_reml, at)
    expect_equal(out_reml, result, tol=1e-6)

})

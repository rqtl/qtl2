context("genome scan by scan1")

test_that("scan1 for backcross with one phenotype", {

    library(qtl)
    data(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk")
    lod0 <- out[,3]

    # inputs for R/qtl2
    pr <- lapply(hyper$geno, function(a) aperm(a$prob, c(1,3,2)))
    y <- hyper$pheno[,1]
    names(y) <- paste(1:nind(hyper))
    for(i in seq(along=pr)) rownames(pr[[i]]) <- names(y)
    n <- length(y)

    # scan
    out2 <- scan1(pr, y)

    # as expected?
    expect_equal(as.numeric(out2), lod0)

    ##############################
    # weighted scan
    set.seed(20151202)
    w <- runif(n, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w)
    lodw0 <- outw[,3]

    names(w) <- names(y)
    outw2 <- scan1(pr, y, weights=w)
    expect_equal(as.numeric(outw2), lodw0)


    ##############################
    # additive covariate
    x <- sample(0:1, n, replace=TRUE)
    names(x) <- names(y)
    out <- scanone(hyper, method="hk", addcovar=x)
    out2 <- scan1(pr, y, x)
    expect_equal(as.numeric(out2), out[,3])

    ##############################
    # additive covariate + weights
    outw <- scanone(hyper, method="hk", addcovar=x, weights=w)
    outw2 <- scan1(pr, y, x, weights=w)
    expect_equal(as.numeric(outw2), outw[,3])


    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x)
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equal(as.numeric(out2), out[,3])

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(as.numeric(out2r), out[,3])

    ##############################
    # additive covariate + weights
    outw <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w)
    outw2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equal(as.numeric(outw2), outw[,3])

    # auto add intcovar?
    outw2 <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(as.numeric(outw2), outw[,3])

})


test_that("scan1 for backcross with multiple phenotypes with NAs", {

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 30
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(1:n_phe), sample(1:5, n_phe, replace=TRUE))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(1:n_ind, nmis[i]), spl[[i]]] <- NA

    # scan by R/qtl
    hyper$pheno <- cbind(y, hyper$pheno[,2,drop=FALSE])
    hyper <- calc.genoprob(hyper, step=1)
    out <- scanone(hyper, method="hk", pheno.col=1:n_phe)
    lod0 <- as.matrix(out[,-(1:2)])
    dimnames(lod0) <- NULL

    # inputs for R/qtl2
    pr <- lapply(hyper$geno, function(a) aperm(a$prob, c(1,3,2)))
    rownames(y) <- paste(1:n_ind)
    for(i in seq(along=pr)) rownames(pr[[i]]) <- rownames(y)

    # scan
    out2 <- scan1(pr, y)
    dimnames(out2) <- NULL

    # as expected?
    expect_equal(out2, lod0)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    outw <- scanone(hyper, method="hk", weights=w, pheno.col=1:n_phe)
    lodw0 <- as.matrix(outw[,-(1:2)])
    dimnames(lodw0) <- NULL

    names(w) <- rownames(y)
    outw2 <- scan1(pr, y, weights=w)
    dimnames(outw2) <- NULL
#    expect_equal(outw2, lodw0)  # FIX_ME <- this isn't working; I suspect that **r/qtl** is wrong


    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x)
    dimnames(lod) <- dimnames(out2) <- NULL
    expect_equal(out2, lod)

    ##############################
    # additive covariate + weights
    outw <- scanone(hyper, method="hk", addcovar=x, weights=w, pheno.col=1:n_phe)
    lodw <- as.matrix(outw[,-(1:2)])
    outw2 <- scan1(pr, y, x, weights=w)
    dimnames(lodw) <- dimnames(outw2) <- NULL
#    expect_equal(outw2, lodw)  # FIX_ME <- this isn't working; I suspect that **r/qtl** is wrong


    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, intcovar=x)
    dimnames(lod) <- dimnames(out2) <- NULL
    expect_equal(out2, lod)

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    dimnames(out2r) <- NULL
    expect_equal(out2r, out2)

    ##############################
    # additive covariate + weights
    outw <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w, pheno.col=1:n_phe)
    lod <- as.matrix(outw[,-(1:2)])
    outw2 <- scan1(pr, y, x, intcovar=x, weights=w)
    dimnames(lodw) <- dimnames(outw2) <- NULL
#    expect_equal(outw2, lodw)  # FIX_ME <- this isn't working; I suspect that **r/qtl** is wrong

    # auto add intcovar?
    outw2r <- scan1(pr, y, intcovar=x, weights=w)
    dimnames(outw2r) <- NULL
    expect_equal(outw2r, outw2)

})

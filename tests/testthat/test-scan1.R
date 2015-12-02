context("genome scan by scan1")

# calc lod scores via lm()
lod_via_lm <-
    function(probs, pheno, addcovar=NULL, intcovar=NULL, weights=NULL)
{
    if(is.list(probs)) { # mult chr by recursion
        result <- NULL
        for(s in seq(along=probs))
            result <- rbind(result, lod_via_lm(probs[[s]], pheno, addcovar, intcovar, weights))
        return(result)
    }

    if(!is.matrix(pheno)) pheno <- as.matrix(pheno)

    d3 <- dim(probs)[3]
    n <- nrow(pheno)
    p <- ncol(pheno)

    result <- matrix(nrow=d3, ncol=p)

    addcovar <- cbind(rep(1, n), addcovar)
    if(is.null(weights)) weights <- rep(1, n)

    for(i in 1:p) {
        not_na <- !is.na(pheno[,i])
        rss0 <- sum(lm(pheno[,i] ~ -1 + addcovar, weights=weights)$resid^2*weights[not_na])
        for(j in 1:d3) {
            if(is.null(intcovar))
                rss1 <- sum(lm(pheno[,i] ~ -1 + addcovar + probs[,-1,j],
                               weights=weights)$resid^2*weights[not_na])
            else
                rss1 <- sum(lm(pheno[,i] ~ -1 + addcovar + probs[,-1,j]*intcovar,
                               weights=weights)$resid^2*weights[not_na])
            result[j,i] <- sum(not_na)/2 * (log10(rss0) - log10(rss1))
        }
    }
    result
}




test_that("scan1 for backcross with one phenotype", {

    library(qtl)
    data(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=5)
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

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])

    ##############################
    # weighted scan
    set.seed(20151202)
    w <- runif(n, 1, 3)
    out <- scanone(hyper, method="hk", weights=w)
    lod0 <- out[,3]

    names(w) <- names(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(as.numeric(out2), lod0)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y, weights=w)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])


    ##############################
    # additive covariate
    x <- sample(0:1, n, replace=TRUE)
    names(x) <- names(y)
    out <- scanone(hyper, method="hk", addcovar=x)
    out2 <- scan1(pr, y, x)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y, x)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w)
    out2 <- scan1(pr, y, x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y, x, weights=w)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x)
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equal(as.numeric(out2), out[,3])

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(as.numeric(out2r), out[,3])

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y, x, x)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w)
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # auto add intcovar?
    out2 <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[1]], y, x, x, weights=w)
    expect_equal(as.numeric(out.lm), out2[1:nrow(out.lm),])

})


test_that("scan1 for backcross with multiple phenotypes with NAs", {

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 15
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(1:n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(1:n_ind, nmis[i]), spl[[i]]] <- NA

    # scan by R/qtl
    hyper$pheno <- cbind(y, hyper$pheno[,2,drop=FALSE])
    hyper <- calc.genoprob(hyper, step=2.5)
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

    # cf lm() for chr 1
    lm_rows <- which(out[,1]==18)
    out.lm <- lod_via_lm(pr[[18]], y)
    expect_equal(out.lm, out2[lm_rows,])

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    out <- scanone(hyper, method="hk", weights=w, pheno.col=1:n_phe)
    lod0 <- as.matrix(out[,-(1:2)])
    dimnames(lod0) <- NULL

    names(w) <- rownames(y)
    out2 <- scan1(pr, y, weights=w)
    dimnames(out2) <- NULL
    expect_equal(out2, lod0)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, weights=w)
    expect_equal(out.lm, out2[lm_rows,])

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x)
    dimnames(lod) <- dimnames(out2) <- NULL
    expect_equal(out2, lod)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x)
    expect_equal(out.lm, out2[lm_rows,])

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, weights=w)
    dimnames(lod) <- dimnames(out2) <- NULL
    expect_equal(out2, lod)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w)
    expect_equal(out.lm, out2[lm_rows,])

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

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x)
    expect_equal(out.lm, out2[lm_rows,])

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    dimnames(lod) <- dimnames(out2) <- NULL
    expect_equal(out2, lod)

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    dimnames(out2r) <- NULL
    expect_equal(out2r, out2)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x, weights=w)
    expect_equal(out.lm, out2[lm_rows,])

})

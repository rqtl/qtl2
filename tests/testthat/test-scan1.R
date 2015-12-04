context("genome scan by scan1")

# calc lod scores via lm(), just one chromosome
lod_via_lm <-
    function(probs, pheno, addcovar=NULL, intcovar=NULL, weights=NULL)
{
    if(!is.matrix(pheno)) pheno <- as.matrix(pheno)

    d3 <- dim(probs)[3]
    n <- nrow(pheno)
    p <- ncol(pheno)

    result <- matrix(nrow=d3, ncol=p)

    addcovar <- cbind(rep(1, n), addcovar)
    if(is.null(weights)) weights <- rep(1, n)

    sample_size <- rep(NA, ncol(pheno))
    names(sample_size) <- colnames(pheno)

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
        sample_size[i] <- sum(not_na)
    }

    attr(result, "sample_size") <- sample_size
    attr(result, "addcovar") <- colnames4attr(addcovar[,-1]) # drop intercept
    attr(result, "intcovar") <- colnames4attr(intcovar)
    dimnames(result) <- list(dimnames(probs)[[3]], colnames(pheno))

    result
}

subset_scan1result <-
    function(x, rows)
{
    at <- attributes(x)
    x <- x[rows,,drop=FALSE]
    for(i in c("addcovar", "intcovar", "sample_size"))
        attr(x, i) <- at[[i]]
    x
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
    expect_equivalent(as.numeric(out2), lod0)

    # cf lm() for chr 18
    lm_rows <- which(out[,1]==18)
    out.lm <- lod_via_lm(pr[[18]], y)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # weighted scan
    set.seed(20151202)
    w <- runif(n, 1, 3)
    out <- scanone(hyper, method="hk", weights=w)
    lod0 <- out[,3]

    names(w) <- names(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(as.numeric(out2), lod0)

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)


    ##############################
    # additive covariate
    x <- sample(0:1, n, replace=TRUE)
    names(x) <- names(y)
    out <- scanone(hyper, method="hk", addcovar=x)
    out2 <- scan1(pr, y, x)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w)
    out2 <- scan1(pr, y, x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x)
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equal(as.numeric(out2), out[,3])

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(as.numeric(out2r), out[,3])

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w)
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # auto add intcovar?
    out2 <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(as.numeric(out2), out[,3])

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

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

    # as expected?
    expect_equivalent(out2, lod0)

    # cf lm() for chr 1
    lm_rows <- which(out[,1]==18)
    out.lm <- lod_via_lm(pr[[18]], y)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    out <- scanone(hyper, method="hk", weights=w, pheno.col=1:n_phe)
    lod0 <- as.matrix(out[,-(1:2)])

    names(w) <- rownames(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equivalent(out2, lod0)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x)
    expect_equivalent(out2, lod)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, weights=w)
    expect_equivalent(out2, lod)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equivalent(out2, lod)

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equivalent(out2r, out2)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w, pheno.col=1:n_phe)
    lod <- as.matrix(out[,-(1:2)])
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equivalent(out2, lod)

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    expect_equivalent(out2r, out2)

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

})

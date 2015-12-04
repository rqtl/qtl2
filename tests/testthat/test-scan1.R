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
    if(is.null(weights)) wts <- rep(1, n)
    else wts <- weights

    sample_size <- rep(NA, ncol(pheno))
    names(sample_size) <- colnames(pheno)

    for(i in 1:p) {
        not_na <- !is.na(pheno[,i])
        rss0 <- sum(lm(pheno[,i] ~ -1 + addcovar, weights=wts)$resid^2*wts[not_na])
        for(j in 1:d3) {
            if(is.null(intcovar))
                rss1 <- sum(lm(pheno[,i] ~ -1 + addcovar + probs[,-1,j],
                               weights=wts)$resid^2*wts[not_na])
            else
                rss1 <- sum(lm(pheno[,i] ~ -1 + addcovar + probs[,-1,j]*intcovar,
                               weights=wts)$resid^2*wts[not_na])
            result[j,i] <- sum(not_na)/2 * (log10(rss0) - log10(rss1))
        }
        sample_size[i] <- sum(not_na)
    }

    attr(result, "sample_size") <- sample_size
    attr(result, "addcovar") <- colnames4attr(addcovar[,-1]) # drop intercept
    attr(result, "intcovar") <- colnames4attr(intcovar)
    if(!is.null(weights)) attr(result, "weights") <- TRUE
    dimnames(result) <- list(dimnames(probs)[[3]], colnames(pheno))

    result
}

# subset rows of scan1 results but preserve attributes
subset_scan1result <-
    function(x, rows)
{
    at <- attributes(x)
    x <- x[rows,,drop=FALSE]
    for(i in c("addcovar", "intcovar", "sample_size", "weights"))
        attr(x, i) <- at[[i]]
    x
}

# revise scanone results to be as expected from scan1
scanone2scan1 <-
    function(x, n=NULL, posnames=NULL, phenames=NULL, addcovar=NULL, intcovar=NULL, Xcovar=NULL, weights=NULL)
{
    x <- as.matrix(x[,-(1:2), drop=FALSE])
    dimnames(x) <- list(posnames, phenames)

    attr(x, "sample_size") <- n
    attr(x, "addcovar") <- colnames4attr(addcovar)
    attr(x, "intcovar") <- colnames4attr(intcovar)
    attr(x, "Xcovar") <- colnames4attr(Xcovar)
    attr(x, "weights") <- weights
    x
}

# now finally to some tests
test_that("scan1 for backcross with one phenotype", {

    library(qtl)
    data(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=5)
    out <- scanone(hyper, method="hk")

    # inputs for R/qtl2
    pr <- lapply(hyper$geno, function(a) aperm(a$prob, c(1,3,2)))
    y <- hyper$pheno[,1]
    names(y) <- paste(1:nind(hyper))
    for(i in seq(along=pr)) rownames(pr[[i]]) <- names(y)
    n <- length(y)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))

    # scan
    out2 <- scan1(pr, y)

    # as expected?
    expect_equal(out2, scanone2scan1(out, 250, posnames))

    # cf lm() for chr 18
    lm_rows <- which(out[,1]==18)
    out.lm <- lod_via_lm(pr[[18]], y)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # weighted scan
    set.seed(20151202)
    w <- runif(n, 1, 3)
    out <- scanone(hyper, method="hk", weights=w)

    names(w) <- names(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, weights=TRUE))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)


    ##############################
    # additive covariate
    x <- sample(0:1, n, replace=TRUE)
    names(x) <- names(y)
    out <- scanone(hyper, method="hk", addcovar=x)
    out2 <- scan1(pr, y, x)
    expect_equal(out2, scanone2scan1(out, 250, posnames, addcovar=x))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w)
    out2 <- scan1(pr, y, x, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, weights=TRUE))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x)
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, intcovar=x))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(out2r, scanone2scan1(out, 250, posnames, colnames(y),
                                      addcovar=x, intcovar=x))


    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w)
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, intcovar=x, weights=TRUE))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(out2r, scanone2scan1(out, 250, posnames, colnames(y),
                                      addcovar=x, intcovar=x, weights=TRUE))

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

    # inputs for R/qtl2
    pr <- lapply(hyper$geno, function(a) aperm(a$prob, c(1,3,2)))
    rownames(y) <- paste(1:n_ind)
    for(i in seq(along=pr)) rownames(pr[[i]]) <- rownames(y)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))

    # scan
    out2 <- scan1(pr, y)

    # as expected?
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)), posnames,
                                     colnames(y)))

    # cf lm() for chr 1
    lm_rows <- which(out[,1]==18)
    out.lm <- lod_via_lm(pr[[18]], y)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    out <- scanone(hyper, method="hk", weights=w, pheno.col=1:n_phe)

    names(w) <- rownames(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)), posnames,
                                     colnames(y), weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, x)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, x, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, x, intcovar=x)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     intcovar=x))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(out2r, scanone2scan1(out, colSums(!is.na(y)),
                                      posnames, colnames(y), addcovar=x,
                                      intcovar=x))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, x, intcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     intcovar=x, weights=TRUE))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(out2r, scanone2scan1(out, colSums(!is.na(y)),
                                      posnames, colnames(y), addcovar=x,
                                      intcovar=x, weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x, weights=w)
    expect_equal(subset_scan1result(out2, lm_rows), out.lm)

})

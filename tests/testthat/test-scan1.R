context("genome scan by scan1")

isx_from_map <-
    function(map)
{
    isx <- rep(FALSE, length(map))
    names(isx) <- names(map)
    isx["X"] <- TRUE

    isx
}


# calc lod scores via lm(), just one chromosome
lod_via_lm <-
    function(probs, pheno, addcovar=NULL, intcovar=NULL, weights=NULL, map)
{
    if(!is.matrix(pheno)) pheno <- as.matrix(pheno)
    if(is.null(colnames(pheno)))
        colnames(pheno) <- paste0("pheno", 1:ncol(pheno))

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
        if(!is.null(addcovar)) not_na <- not_na & complete.cases(addcovar)
        if(!is.null(intcovar)) not_na <- not_na & complete.cases(intcovar)

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

    dimnames(result) <- list(dimnames(probs)[[3]], colnames(pheno))

    attr(result, "sample_size") <- sample_size

    class(result) <- c("scan1", "matrix")
    result
}

# revise scanone results to be as expected from scan1
scanone2scan1 <-
    function(x, n=NULL, posnames=NULL, phenames=NULL, addcovar=NULL, intcovar=NULL, Xcovar=NULL, weights=FALSE)
{
    map <- split(x[,2], factor(x[,1], levels=unique(x[,1])))
    names(map) <- unique(x[,1])
    if(!is.null(posnames)) {
        mapn <- split(posnames, factor(x[,1], levels=unique(x[,1])))
        for(i in seq(along=map)) names(map[[i]]) <- mapn[[i]]
    }

    x <- as.matrix(x[,-(1:2), drop=FALSE])
    if(is.null(phenames)) phenames <- paste0("pheno", 1:ncol(x))
    dimnames(x) <- list(posnames, phenames)

    if(is.null(names(n))) names(n) <- phenames

    result <- x
    attr(result, "sample_size") <- n

    class(result) <- c("scan1", "matrix")
    result
}

convert_probs2qtl2 <-
    function(cross)
{
    result <- lapply(cross$geno, function(a) {
        pr <- aperm(a$prob, c(1,3,2))
        rownames(pr) <- paste(1:nrow(pr))
        pr })

    attr(result, "is_x_chr") <- isx_from_map(convert2cross2(cross)$gmap)
    attr(result, "crosstype") <- class(cross)[1]

    class(result) <- c("calc_genoprob", "list")
    result
}

# now finally to some tests
test_that("scan1 for backcross with one phenotype", {

    library(qtl)
    data(hyper)

    # scan by R/qtl
    hyper <- calc.genoprob(hyper, step=5)
    out <- scanone(hyper, method="hk")

    # inputs for R/qtl2
    pr <- convert_probs2qtl2(hyper)
    y <- hyper$pheno[,1]
    names(y) <- paste(1:nind(hyper))
    for(i in seq(along=pr)) rownames(pr[[i]]) <- names(y)
    n <- length(y)

    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    # scan
    out2 <- scan1(pr, y)

    # as expected?
    expect_equal(out2, scanone2scan1(out, 250, posnames, weights=FALSE))

    map <- lapply(hyper$geno, function(a) attr(a$prob, "map"))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # weighted scan
    set.seed(20151202)
    w <- runif(n, 1, 3)
    out <- scanone(hyper, method="hk", weights=w)

    names(w) <- names(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, weights=TRUE))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)


    ##############################
    # additive covariate
    x <- sample(0:1, n, replace=TRUE)
    names(x) <- names(y)
    out <- scanone(hyper, method="hk", addcovar=x)
    out2 <- scan1(pr, y, addcovar=x)
    expect_equal(out2, scanone2scan1(out, 250, posnames, addcovar=x))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w)
    out2 <- scan1(pr, y, addcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, weights=TRUE))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x)
    out2 <- scan1(pr, y, addcovar=x, intcovar=x)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, intcovar=x))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(out2r, scanone2scan1(out, 250, posnames, colnames(y),
                                      addcovar=x, intcovar=x))


    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, x, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w)
    out2 <- scan1(pr, y, addcovar=x, intcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, 250, posnames, colnames(y),
                                     addcovar=x, intcovar=x, weights=TRUE))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(out2r, scanone2scan1(out, 250, posnames, colnames(y),
                                      addcovar=x, intcovar=x, weights=TRUE))

    # cf lm() for chr 18
    out.lm <- lod_via_lm(pr[[18]], y, x, x, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)

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
    spl <- split(sample(n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(n_ind, nmis[i]), spl[[i]]] <- NA

    # scan by R/qtl
    hyper$pheno <- cbind(y, hyper$pheno[,2,drop=FALSE])
    hyper <- calc.genoprob(hyper, step=2.5)
    out <- scanone(hyper, method="hk", pheno.col=1:n_phe)

    map <- lapply(hyper$geno, function(a) attr(a$prob, "map"))

    # inputs for R/qtl2
    pr <- convert_probs2qtl2(hyper)
    rownames(y) <- paste(1:n_ind)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    # scan
    out2 <- scan1(pr, y)

    # as expected?
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)), posnames,
                                     colnames(y)))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    out <- scanone(hyper, method="hk", weights=w, pheno.col=1:n_phe)

    names(w) <- rownames(y)
    out2 <- scan1(pr, y, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)), posnames,
                                     colnames(y), weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, addcovar=x)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # additive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, weights=w, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, addcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # interactive covariate
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, addcovar=x, intcovar=x)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     intcovar=x))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x)
    expect_equal(out2r, scanone2scan1(out, colSums(!is.na(y)),
                                      posnames, colnames(y), addcovar=x,
                                      intcovar=x))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x, map=map)
    expect_equal(out2[map,"18",], out.lm)

    ##############################
    # interactive covariate + weights
    out <- scanone(hyper, method="hk", addcovar=x, intcovar=x, weights=w, pheno.col=1:n_phe)
    out2 <- scan1(pr, y, addcovar=x, intcovar=x, weights=w)
    expect_equal(out2, scanone2scan1(out, colSums(!is.na(y)),
                                     posnames, colnames(y), addcovar=x,
                                     intcovar=x, weights=TRUE))

    # auto add intcovar?
    out2r <- scan1(pr, y, intcovar=x, weights=w)
    expect_equal(out2r, scanone2scan1(out, colSums(!is.na(y)),
                                      posnames, colnames(y), addcovar=x,
                                      intcovar=x, weights=TRUE))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, intcovar=x, weights=w, map=map)
    expect_equal(out2[map,"18",], out.lm)

})


test_that("scan1 works with NAs in the covariates", {

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 15
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(n_ind, nmis[i]), spl[[i]]] <- NA
    hyper$pheno <- cbind(y, hyper$pheno[,2,drop=FALSE])

    # genoprobs by R/qtl
    hyper <- calc.genoprob(hyper, step=2.5)

    # inputs for R/qtl2
    pr <- convert_probs2qtl2(hyper)
    rownames(y) <- paste(1:n_ind)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    x[5] <- NA

    suppressWarnings(out <- scanone(hyper, method="hk", addcovar=x, pheno.col=1:n_phe))
    out2 <- scan1(pr, y, addcovar=x)
    expect_equal(out2, scanone2scan1(out, colSums(!(is.na(y) | is.na(x))),
                                     posnames, colnames(y), addcovar=x))

    # cf lm() for chr 1
    out.lm <- lod_via_lm(pr[[18]], y, x, map=map)
    expect_equal(out2[map,"18",], out.lm)
})


test_that("scan1 aligns the individuals", {

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 15
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(n_ind, nmis[i]), spl[[i]]] <- NA
    hyper$pheno <- cbind(y, hyper$pheno[,2,drop=FALSE])

    # genoprobs from R/qtl
    hyper <- calc.genoprob(hyper, step=2.5)

    # inputs for R/qtl2
    pr <- convert_probs2qtl2(hyper)
    rownames(y) <- paste(1:n_ind)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    # scan
    out <- scan1(pr, y)
    out_perm <- scan1(pr, y[sample(n_ind),])
    expect_equal(out_perm, out)

    class(pr) <- c("calc_genoprob", "list") # allows simpler reordering of individuals
    out_perm <- scan1(pr[sample(n_ind),], y)
    expect_equal(out_perm, out)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    names(w) <- rownames(y)

    out <- scan1(pr, y, weights=w)
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], weights=w[sample(n_ind)])
    expect_equal(out_perm, out)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)

    out <- scan1(pr, y, addcovar=x)
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], addcovar=x[sample(n_ind)])
    expect_equal(out_perm, out)

    ##############################
    # additive covariate + weights
    out <- scan1(pr, y, addcovar=x, weights=w)
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], addcovar=x[sample(n_ind)],
                      weights=w[sample(n_ind)])
    expect_equal(out_perm, out)

    ##############################
    # interactive covariate
    out <- scan1(pr, y, addcovar=x, intcovar=x)
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], addcovar=x[sample(n_ind)],
                      intcovar=x[sample(n_ind)])
    expect_equal(out_perm, out)

    # auto add intcovar?
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], intcovar=x[sample(n_ind)])
    expect_equal(out_perm, out)

    ##############################
    # interactive covariate + weights
    out <- scan1(pr, y, addcovar=x, intcovar=x, weights=w)
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),], addcovar=x[sample(n_ind)],
                      intcovar=x[sample(n_ind)], weights=w[sample(n_ind)])
    expect_equal(out_perm, out)

    # auto add intcovar?
    out_perm <- scan1(pr[sample(n_ind),], y[sample(n_ind),],
                      intcovar=x[sample(n_ind)], weights=w[sample(n_ind)])
    expect_equal(out_perm, out)

})


test_that("multi-core scan1 works", {
    if(isnt_karl()) skip("this test only run locally")

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 15
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(n_ind, nmis[i]), spl[[i]]] <- NA
    rownames(y) <- paste(1:n_ind)

    # inputs for R/qtl2
    library(qtl2geno)
    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2, step=2.5)
    pr <- calc_genoprob(hyper2, map)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    # scan
    out <- scan1(pr, y)
    out_multicore <- scan1(pr, y, cores=4)
    expect_equal(out_multicore, out)
    out_multicore <- scan1(pr, y, cores=0) # maximum cores
    expect_equal(out_multicore, out)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    names(w) <- rownames(y)

    out <- scan1(pr, y, weights=w)
    out_multicore <- scan1(pr, y, weights=w, cores=4)
    expect_equal(out_multicore, out)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)

    out <- scan1(pr, y, addcovar=x)
    out_multicore <- scan1(pr, y, addcovar=x, cores=4)
    expect_equal(out_multicore, out)

    ##############################
    # additive covariate + weights
    out <- scan1(pr, y, addcovar=x, weights=w)
    out_multicore <- scan1(pr, y, addcovar=x, weights=w, cores=4)
    expect_equal(out_multicore, out)

    ##############################
    # interactive covariate
    out <- scan1(pr, y, addcovar=x, intcovar=x)
    out_multicore <- scan1(pr, y, addcovar=x, intcovar=x, cores=4)
    expect_equal(out_multicore, out)

    # auto add intcovar?
    out_multicore <- scan1(pr, y, intcovar=x, cores=4)
    expect_equal(out_multicore, out)

    ##############################
    # interactive covariate + weights
    out <- scan1(pr, y, addcovar=x, intcovar=x, weights=w)
    out_multicore <- scan1(pr, y, addcovar=x, intcovar=x, weights=w, cores=4)
    expect_equal(out_multicore, out)

    # auto add intcovar?
    out_multicore <- scan1(pr, y, intcovar=x, weights=w, cores=4)
    expect_equal(out_multicore, out)

})


test_that("scan1 LOD results don't depend on scale of x and y", {

    set.seed(20151202)
    library(qtl)
    data(hyper)

    # phenotypes
    n_phe <- 15
    n_ind <- nind(hyper)
    y <- matrix(rnorm(n_ind*n_phe), ncol=n_phe)
    # 5 batches
    spl <- split(sample(n_phe), rep(1:5, 3))
    nmis <- c(0, 5, 10, 15, 20)
    for(i in seq(along=spl)[-1])
        y[sample(n_ind, nmis[i]), spl[[i]]] <- NA
    rownames(y) <- paste(1:n_ind)

    # inputs for R/qtl2
    library(qtl2geno)
    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=2.5)
    pr <- calc_genoprob(hyper2, map)
    posnames <- unlist(lapply(pr, function(a) dimnames(a)[[3]]))
    names(posnames) <- NULL

    ybig <- y*100

    # scan
    out <- scan1(pr, y)
    outbig <- scan1(pr, ybig)
    expect_equal(outbig, out)

    ##############################
    # weighted scan
    w <- runif(n_ind, 1, 3)
    names(w) <- rownames(y)
    wbig <- w*5

    out <- scan1(pr, y, weights=w)
    outbig <- scan1(pr, ybig, weights=wbig)
    expect_equal(outbig, out)

    ##############################
    # additive covariate
    x <- sample(0:1, n_ind, replace=TRUE)
    names(x) <- rownames(y)
    xbig <- x*50

    out <- scan1(pr, y, addcovar=x)
    outbig <- scan1(pr, ybig, addcovar=xbig)
    expect_equal(outbig, out)

    ##############################
    # additive covariate + weights
    out <- scan1(pr, y, addcovar=x, weights=w)
    outbig <- scan1(pr, ybig, addcovar=xbig, weights=wbig)
    expect_equal(outbig, out)

    ##############################
    # interactive covariate
    out <- scan1(pr, y, addcovar=x, intcovar=x)
    outbig <- scan1(pr, ybig, addcovar=xbig, intcovar=xbig)
    expect_equal(outbig, out)

    ##############################
    # interactive covariate + weights
    out <- scan1(pr, y, addcovar=x, intcovar=x, weights=w)
    outbig <- scan1(pr, ybig, addcovar=xbig, intcovar=xbig, weights=wbig)
    expect_equal(outbig, out)

})

test_that("scan1 deals with mismatching individuals", {
    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs, "loco")
    Xc <- get_x_covar(iron)
    X <- match(iron$covar$sex, c("f", "m"))-1
    names(X) <- rownames(iron$covar)

    ind <- c(1:50, 101:150)
    expected <- scan1(probs[ind,], iron$pheno[ind,,drop=FALSE], addcovar=X[ind], intcovar=X[ind],
                      Xcovar=Xc[ind,])
    expect_equal(scan1(probs[ind,], iron$pheno, addcovar=X, intcovar=X, Xcovar=Xc), expected)
    expect_equal(scan1(probs, iron$pheno[ind,], addcovar=X, intcovar=X, Xcovar=Xc), expected)
    expect_equal(scan1(probs, iron$pheno, addcovar=X[ind], intcovar=X, Xcovar=Xc), expected)
    expect_equal(scan1(probs, iron$pheno, addcovar=X, intcovar=X[ind], Xcovar=Xc), expected)
    expect_equal(scan1(probs, iron$pheno, addcovar=X, intcovar=X, Xcovar=Xc[ind,]), expected)

})

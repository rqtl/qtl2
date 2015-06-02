context("linear regression")
library(qtl)

test_that("lin regr works for simple example", {

    set.seed(27343534)

    n <- 1000
    x <- rnorm(n, 50, 10)
    X <- cbind(1, x)
    y <- 30 + 0.5*x + rnorm(n, 0, 2.5)
    Y <- as.matrix(y)

    lm.out <- lm(y ~ x)
    resid <- lm.out$resid
    names(resid) <- NULL
    rss <- sum(resid^2)

    # cholesky-based regression
    lltFit <- fit_linreg_eigenchol(X, y)
    expect_equal(rss, lltFit$rss)
    expect_equivalent(lm.out$coef, lltFit$coef)
    expect_equivalent(lm.out$fitted, lltFit$fitted)
    expect_equivalent(summary(lm.out)$coef[,2], lltFit$SE)

    lltRSS <- calc_rss_eigenchol(X, y)
    expect_equal(rss, lltRSS)

    # QR-based regression
    qrFit <- fit_linreg_eigenqr(X, y)
    expect_equal(rss, qrFit$rss)
    expect_equivalent(lm.out$coef, qrFit$coef)
    expect_equivalent(lm.out$fitted, qrFit$fitted)
    expect_equivalent(summary(lm.out)$coef[,2], qrFit$SE)

    qrRSS <- calc_rss_eigenqr(X, y)
    expect_equal(rss, qrRSS)

    # LAPACK rss
    lapack_rss <- calc_rss_lapack(X,Y)
    expect_equal(rss, lapack_rss)

    dgelsy_rss <- calc_rss_lapack(X,Y, skip_dgels=TRUE)
    expect_equal(rss, dgelsy_rss)

    # LAPACK resid
    lapack_resid <- calc_resid_lapack(X,Y)
    expect_equal(lapack_resid, as.matrix(resid))

    dgelsy_resid <- calc_resid_lapack(X,Y, skip_dgels=TRUE)
    expect_equal(dgelsy_resid, as.matrix(resid))

})

test_that("lin regr works for reduced-rank example", {

    dd <- data.frame(f1 = gl(4, 6, labels = LETTERS[1:4]),
                     f2 = gl(3, 2, labels = letters[1:3]))[-(7:8), ]
    mm <- model.matrix(~ f1*f2, dd)
    y <- mm %*% seq_len(ncol(mm)) + rnorm(nrow(mm), sd = 0.1)
    Y <- as.matrix(y)

    lm.out <- lm(y ~ -1 + mm)
    resid <- lm.out$resid
    names(resid) <- NULL
    rss <- sum(resid^2)
    lm_se <- lm.out$coef
    lm_se[!is.na(lm_se)] <- summary(lm.out)$coef[,2]

    # (cholesky-based regression won't work here)

    # QR-based regression
    qrFit <- fit_linreg_eigenqr(mm, y)
    expect_equal(rss, qrFit$rss)
    expect_equivalent(lm.out$coef, qrFit$coef)
    expect_equivalent(lm.out$fitted, qrFit$fitted)
    expect_equivalent(lm_se, qrFit$SE)

    qrRSS <- calc_rss_eigenqr(mm, y)
    expect_equal(rss, qrRSS)

    # LAPACK RSS
    lapack_rss <- calc_rss_lapack(mm,Y)
    expect_equal(rss, lapack_rss)

    dgelsy_rss <- calc_rss_lapack(mm, Y, skip_dgels=TRUE)
    expect_equal(rss, dgelsy_rss)

    # LAPACK resid
    lapack_resid <- calc_resid_lapack(mm,Y)
    expect_equal(lapack_resid, as.matrix(resid))

    dgelsy_resid <- calc_resid_lapack(mm, Y, skip_dgels=TRUE)
    expect_equal(dgelsy_resid, as.matrix(resid))

})


test_that("lin regr works in a QTL situation", {

    data(hyper)
    hyper <- calc.genoprob(hyper, step=1, err=0.001)
    qtl <- makeqtl(hyper, c(1,4,6,15), c(68.3, 30, 66.7, 17.5),
                   what="prob")

    n <- nind(hyper)
    y <- hyper$pheno[,1]
    Y <- as.matrix(y)
    X <- matrix(ncol=6, nrow=n)
    X[,1] <- 1
    for(i in 2:5) X[,i] <- qtl$prob[[i-1]][,2]
    X[,6] <- X[,4]*X[,5]

    lodfull <- fitqtl(hyper, qtl=qtl, method="hk", dropone=FALSE,
                      formula=y~q1+q2+q3*q4)$lod
    lodadd <- fitqtl(hyper, qtl=qtl, method="hk", dropone=FALSE,
                     formula=y~q1+q2+q3+q4)$lod
    lod14 <- fitqtl(hyper, qtl=qtl, method="hk", dropone=FALSE,
                    formula=y~q1+q2)$lod
    lod4 <- fitqtl(hyper, qtl=qtl, method="hk", dropone=FALSE,
                    formula=y~q2)$lod

    # eigen qr
    rss0 <- calc_rss_eigenqr(X[,1,drop=FALSE], y)
    expect_equal(rss0, sum(lm(y~1)$resid^2))

    rssfull <- calc_rss_eigenqr(X, y)
    expect_equal(rssfull, sum(lm(y ~ X)$resid^2))
    expect_equal(lodfull, (n/2)*log10(rss0/rssfull))

    rssadd <- calc_rss_eigenqr(X[,1:5], y)
    expect_equal(rssadd, sum(lm(y ~ X[,2:5])$resid^2))
    expect_equal(lodadd, (n/2)*log10(rss0/rssadd))

    rss14 <- calc_rss_eigenqr(X[,1:3], y)
    expect_equal(rss14, sum(lm(y ~ X[,2:3])$resid^2))
    expect_equal(lod14, (n/2)*log10(rss0/rss14))

    rss4 <- calc_rss_eigenqr(X[,c(1,3)], y)
    expect_equal(rss4, sum(lm(y ~ X[,3,drop=FALSE])$resid^2))
    expect_equal(lod4, (n/2)*log10(rss0/rss4))

})

test_that("lin regr works for multiple columns", {

    set.seed(27343534)

    n <- 1000
    x <- rnorm(n, 50, 10)
    X <- cbind(1, x)
    y <- 30 + 0.5*x + rnorm(n, 0, 2.5)

    ncolY <- 10
    Y <- permute_nvector(10, y)

    resid <- lm.fit(X,Y)$resid
    lm.rss <- colSums(resid^2)

    # RSS
    expect_equal(lm.rss, calc_mvrss_eigenchol(X, Y))
    expect_equal(lm.rss, calc_mvrss_eigenqr(X, Y))
    expect_equal(lm.rss, calc_rss_lapack(X, Y))
    expect_equal(lm.rss, calc_rss_lapack(X, Y, skip_dgels=TRUE))

    # residuals
    expect_equal(resid, calc_resid_lapack(X, Y))
    expect_equal(resid, calc_resid_lapack(X, Y, skip_dgels=TRUE))

})

test_that("lin regr works for multiple columns, reduced-rank X", {

    dd <- data.frame(f1 = gl(4, 6, labels = LETTERS[1:4]),
                     f2 = gl(3, 2, labels = letters[1:3]))[-(7:8), ]
    mm <- model.matrix(~ f1*f2, dd)
    y <- mm %*% seq_len(ncol(mm)) + rnorm(nrow(mm), sd = 0.1)
    Y <- permute_nvector(10, y)

    resid <- lm.fit(mm, Y)$resid
    lm.rss <- colSums(resid^2)

    # RSS
    expect_equal(lm.rss, calc_mvrss_eigenqr(mm, Y))
    expect_equal(lm.rss, calc_rss_lapack(mm, Y))
    expect_equal(lm.rss, calc_rss_lapack(mm, Y, skip_dgels=TRUE))

    # resid
    expect_equal(resid, calc_resid_lapack(mm, Y))
    expect_equal(resid, calc_resid_lapack(mm, Y, skip_dgels=TRUE))

})

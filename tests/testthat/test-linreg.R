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

    # QR-based regression
    qrFit <- fit_linreg_eigenqr(X, y, TRUE)
    expect_equal(rss, qrFit$rss)
    expect_equivalent(lm.out$coef, qrFit$coef)
    expect_equivalent(lm.out$fitted, qrFit$fitted)
    expect_equivalent(summary(lm.out)$coef[,2], qrFit$SE)

    qrRSS <- calc_rss_eigenqr(X, y)
    expect_equal(rss, qrRSS)

    # eigen residuals
    eigenqr_resid <- calc_resid_eigenqr(X,Y)
    expect_equal(eigenqr_resid, as.matrix(resid))
    eigenchol_resid <- calc_resid_eigenqr(X,Y)
    expect_equal(eigenchol_resid, as.matrix(resid))

    # generic
    linreg_rss <- calc_rss_linreg(X,Y)
    expect_equal(linreg_rss, rss)
    linreg_resid <- calc_resid_linreg(X,Y)
    expect_equal(linreg_resid, as.matrix(resid))

    # just the coefficients
    coef <- calc_coef_linreg(X,Y)
    expect_equal(coef, as.numeric(lm.out$coef))

    # coefficients and SEs
    coefSE <- calc_coefSE_linreg(X,Y)
    expected <- summary(lm.out)$coef
    expected <- list(coef=as.numeric(expected[,1]),
                     SE=as.numeric(expected[,2]))
    expect_equal(coefSE, expected)

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
    qrFit <- fit_linreg_eigenqr(mm, y, TRUE)
    expect_equal(rss, qrFit$rss)
    expect_equivalent(lm.out$coef, qrFit$coef)
    expect_equivalent(lm.out$fitted, qrFit$fitted)
    expect_equivalent(lm_se, qrFit$SE)

    qrRSS <- calc_rss_eigenqr(mm, y)
    expect_equal(rss, qrRSS)

    # eigen resid
    eigenqr_resid <- calc_resid_eigenqr(mm, Y)
    expect_equal(eigenqr_resid, as.matrix(resid))

    # generic
    linreg_rss <- calc_rss_linreg(mm,Y)
    expect_equal(linreg_rss, rss)
    linreg_resid <- calc_resid_linreg(mm,Y)
    expect_equal(linreg_resid, as.matrix(resid))

    # just the coefficients
    coef <- calc_coef_linreg(mm,y)
    expect_equal(coef, as.numeric(lm.out$coef))

    # coefficients and SEs
    coefSE <- calc_coefSE_linreg(mm,y)
    ### a bit of effort due to one coefficient being NA
    expected <- list(coef=as.numeric(lm.out$coef),
                     SE=as.numeric(lm.out$coef))
    lm.coef <- summary(lm.out)$coef
    for(i in 1:2) expected[[i]][!is.na(coef)] <- lm.coef[,i]
    expect_equal(coefSE, expected)

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

    # residuals
    expect_equal(resid, calc_resid_eigenqr(X, Y))
    expect_equal(resid, calc_resid_eigenchol(X, Y))

    # generic
    expect_equal(lm.rss, calc_rss_linreg(X, Y))
    expect_equal(as.matrix(resid), calc_resid_linreg(X, Y))

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

    # resid
    expect_equal(resid, calc_resid_eigenqr(mm, Y))

    # generic
    expect_equal(lm.rss, calc_rss_linreg(mm, Y))
    expect_equal(as.matrix(resid), calc_resid_linreg(mm, Y))
})


test_that("calculation of residuals for 3d arrays works", {

    library(qtl)
    data(hyper)
    hyper <- hyper[1,]
    hyper2 <- qtl2geno::convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1)
    pr <- qtl2geno::calc_genoprob(hyper2, map, error_prob=0.002)
    pr <- aperm(pr[[1]], c(1,3,2)) # reorient to have genomic position last

    # residuals with intercept plus the phenotype
    X <- cbind(1, hyper$pheno[,1])
    expected <- array(0, dim=dim(pr))
    for(i in 1:dim(pr)[3])
        expected[,,i] <- lm(pr[,,i] ~ X)$resid

    resid <- calc_resid_linreg_3d(X, pr)
    expect_equal(expected, resid)

    # residuals with just the intercept
    expected <- array(0, dim=dim(pr))
    for(i in 1:dim(pr)[3])
        expected[,,i] <- lm(pr[,,i] ~ 1)$resid

    resid <- calc_resid_linreg_3d(X[,1,drop=FALSE], pr)
    expect_equal(expected, resid)

})

context("fit1")

test_that("fit1 by H-K works in intercross", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(18:19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov) })

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(3, 4, 12)
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov) })

    pos <- cumsum(c(0, npos[-3])) + pmar
    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:3) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:3)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # repeat the whole thing with a couple of missing phenotypes
    pheno[c(187, 244)] <- NA

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, se=TRUE) })

    # fit1, no missing data
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, se=TRUE) })

    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:3) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:3)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # check SEs
    for(i in 1:3)
        expect_equal(out_fit1[[i]]$SE, attr(coef[[i]], "SE")[pmar[i],])

    # direct calculations, chr 18
    lm0 <- lm(pheno ~ covar)
    X <- cbind(probs[[1]][,,pmar[1]], covar)
    colnames(X) <- c("SS", "SB", "BB", "ac1")
    lm1 <- lm(pheno ~ -1 + X)
    n <- sum(!is.na(pheno))
    sigsq0 <- sum(lm0$resid^2)/n
    sigsq1 <- sum(lm1$resid^2)/n
    expect_equal(out_fit1[[1]]$lod, n/2*log10(sum(lm0$resid^2)/sum(lm1$resid^2)))
    expect_equal(out_fit1[[1]]$ind_lod, (dnorm(lm1$resid,0,sqrt(sigsq1), TRUE) - dnorm(lm0$resid,0,sqrt(sigsq0),TRUE))/log(10))
    expect_equal(out_fit1[[1]]$coef, stats::setNames(lm1$coef, c("SS", "SB", "BB", "ac1")))
    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(lm1)$coef[,2], c("SS", "SB", "BB", "ac1")))

    # direct calculations, chr X
    lm0 <- lm(pheno ~ Xcovar)
    X <- probs[[3]][,,pmar[3]]
    colnames(X) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    lm1 <- lm(pheno ~ -1 + X)
    n <- sum(!is.na(pheno))
    sigsq0 <- sum(lm0$resid^2)/n
    sigsq1 <- sum(lm1$resid^2)/n
    expect_equal(out_fit1[[3]]$lod, n/2*log10(sum(lm0$resid^2)/sum(lm1$resid^2)))
    expect_equal(out_fit1[[3]]$ind_lod, (dnorm(lm1$resid,0,sqrt(sigsq1), TRUE) - dnorm(lm0$resid,0,sqrt(sigsq0),TRUE))/log(10))
    expect_equal(out_fit1[[3]]$coef, stats::setNames(lm1$coef, c("SS", "SB", "BS", "BB", "SY", "BY")))
    expect_equal(out_fit1[[3]]$SE, stats::setNames(summary(lm1)$coef[,2], c("SS", "SB", "BS", "BB", "SY", "BY")))

})


test_that("fit1 by H-K works in riself", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    grav2 <- grav2[,4:5]
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)

    pheno <- grav2$pheno[,219]

    # calculate LOD scores
    out <- scan1(probs, pheno)

    # estimate coefficients
    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno))

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(1, 172)
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno))

    pos <- c(0,npos[1]) + pmar
    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:2) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:2)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # repeat the whole thing with a couple of missing phenotypes
    pheno[c(24, 106)] <- NA

    # calculate LOD scores
    out <- scan1(probs, pheno)

    # estimate coefficients
    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno, se=TRUE))

    # fit1, no missing data
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno, se=TRUE))

    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:2) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:2)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])


    # check SEs
    for(i in 1:2)
        expect_equal(out_fit1[[i]]$SE, attr(coef[[i]], "SE")[pmar[i],])

    # direct calculations, chr 18
    lm0 <- lm(pheno ~ 1)
    X <- probs[[1]][,,pmar[1]]
    colnames(X) <- c("LL", "CC")
    lm1 <- lm(pheno ~ -1 + X)
    n <- sum(!is.na(pheno))
    sigsq0 <- sum(lm0$resid^2)/n
    sigsq1 <- sum(lm1$resid^2)/n
    expect_equal(out_fit1[[1]]$lod, n/2*log10(sum(lm0$resid^2)/sum(lm1$resid^2)))
    expect_equal(out_fit1[[1]]$ind_lod, (dnorm(lm1$resid,0,sqrt(sigsq1), TRUE) - dnorm(lm0$resid,0,sqrt(sigsq0),TRUE))/log(10))
    expect_equal(out_fit1[[1]]$coef, stats::setNames(lm1$coef, c("LL", "CC")))
    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(lm1)$coef[,2], c("LL", "CC")))

})


test_that("fit1 by LMM works in intercross", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    kinship_loco <- calc_kinship(probs, "loco")
    probs <- probs[,c(7,19,"X")]
    kinship_loco <- kinship_loco[c(7,19,"X")]

    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # calculate LOD scores
    out <- scan1(probs, pheno, kinship, addcovar=covar, Xcovar=Xcovar)
    out_loco <- scan1(probs, pheno, kinship_loco, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients (with SEs)
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship, addcovar=cov, se=TRUE) })

    coef_loco <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship_loco[[i]], addcovar=cov, se=TRUE) })

    names(coef) <- names(coef_loco) <- names(kinship_loco)

    # positions to test fit1()
    chr <- rep(c(7,19,"X"), each=4)
    n_pos <- sapply(probs, function(a) dim(a)[3])
    index <- c(1,6,7,58, 1,15,31,36,  1,20,27,30)
    pos <- index + rep(cumsum(c(0,n_pos)), c(4,4,4,0))

    # fit1, no missing data
    out_fit1 <- lapply(seq(along=chr),
                       function(i) {
        if(chr[i]=="X") { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[chr[i]]][,,index[i]], pheno, kinship, addcovar=cov, nullcovar=nullcov, se=TRUE) })

    out_fit1_loco <- lapply(seq(along=chr),
                       function(i) {
        if(chr[i]=="X") { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[chr[i]]][,,index[i]], pheno, kinship_loco[[chr[i]]], addcovar=cov, nullcovar=nullcov, se=TRUE) })

    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in seq(along=out_fit1)) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1], tol=1e-6)
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }
    # same with _loco version
    for(i in seq(along=out_fit1)) {
        expect_equal(out_fit1_loco[[i]]$lod, out_loco[pos[i],1], tol=1e-6)
        expect_equal(sum(out_fit1_loco[[i]]$ind_lod), out_fit1_loco[[i]]$lod)
    }

    # check coefficients
    for(i in seq(along=out_fit1))
        expect_equal(out_fit1[[i]]$coef, coef[[chr[i]]][index[i],])
    # same, _loco version
    for(i in seq(along=out_fit1))
        expect_equal(out_fit1_loco[[i]]$coef, coef_loco[[chr[i]]][index[i],])

    # check SEs
    for(i in seq(along=out_fit1))
        expect_equal(out_fit1[[i]]$SE, attr(coef[[chr[i]]], "SE")[index[i],])
    # same, _loco version
    for(i in seq(along=out_fit1))
        expect_equal(out_fit1_loco[[i]]$SE, attr(coef_loco[[chr[i]]], "SE")[index[i],])

})

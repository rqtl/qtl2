context("fit1")

test_that("fit1 by H-K works in intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
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

    # fit1, missing data
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


test_that("fit1 by H-K works in intercross, with weights", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(18:19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    w <- setNames(runif(n_ind(iron), 1, 10), ind_ids(iron))

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar, weights=w)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, weights=w) })

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(3, 4, 12)
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, weights=w) })

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
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar, weights=w)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, se=TRUE, weights=w) })

    # fit1, missing data
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, se=TRUE, weights=w) })

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

})


test_that("fit1 by H-K works in riself", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
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

    # fit1, missing data
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

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
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
        if(i==3) { cov <- NULL; nullcov <- Xcovar }
        else { cov <- covar; nullcov <- NULL }
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship,
                  addcovar=cov, nullcovar=nullcov, se=TRUE) })

    coef_loco <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) { cov <- NULL; nullcov <- Xcovar }
        else { cov <- covar; nullcov <- NULL }
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship_loco[i],
                  addcovar=cov, nullcovar=nullcov, se=TRUE) })

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


test_that("fit1 by LMM works in intercross, with weights", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    kinship <- calc_kinship(probs)
    kinship_loco <- calc_kinship(probs, "loco")
    probs <- probs[,c(7,19,"X")]
    kinship_loco <- kinship_loco[c(7,19,"X")]

    set.seed(2661343)
    w <- setNames(runif(n_ind(iron), 1, 1), ind_ids(iron))

    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # calculate LOD scores
    out <- scan1(probs, pheno, kinship, addcovar=covar, Xcovar=Xcovar, weights=w)
    out_loco <- scan1(probs, pheno, kinship_loco, addcovar=covar, Xcovar=Xcovar, weights=w)

    # estimate coefficients (with SEs)
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) { cov <- NULL; nullcov <- Xcovar }
        else { cov <- covar; nullcov <- NULL }
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship,
                  addcovar=cov, nullcovar=nullcov, se=TRUE, weights=w) })

    coef_loco <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) { cov <- NULL; nullcov <- Xcovar }
        else { cov <- covar; nullcov <- NULL }
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, kinship_loco[i],
                  addcovar=cov, nullcovar=nullcov, se=TRUE, weights=w) })

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
        fit1(probs[[chr[i]]][,,index[i]], pheno, kinship, addcovar=cov, nullcovar=nullcov, se=TRUE, weights=w) })

    out_fit1_loco <- lapply(seq(along=chr),
                       function(i) {
        if(chr[i]=="X") { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[chr[i]]][,,index[i]], pheno, kinship_loco[[chr[i]]], addcovar=cov, nullcovar=nullcov,
             se=TRUE, weights=w) })

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




test_that("fit1 handles contrasts properly in a backcross", {

    library(qtl)
    data(fake.bc)
    fake.bc <- subset(fake.bc, chr=c(2,5))

    fake_bc <- convert2cross2(fake.bc)

    fake.bc <- calc.genoprob(fake.bc, error.prob=0.002, map.function="c-f")
    pr <- calc_genoprob(fake_bc, error_prob=0.002, map_function="c-f")

    # fit single-QTL model at two positions
    qtl <- makeqtl(fake.bc, c(2,5), pos=c(44.8, 9.8), what="prob")
    out.fq1 <- fitqtl(fake.bc, qtl=qtl, formula=y~Q1, get.ests=TRUE, method="hk")
    co1 <- summary(out.fq1)$ests
    out.fq2 <- fitqtl(fake.bc, qtl=qtl, formula=y~Q2, get.ests=TRUE, method="hk")
    co2 <- summary(out.fq2)$ests

    contrasts <- cbind(mu=c(1,1), a=c(-0.5,0.5))
    out_fit1a <- fit1(pull_genoprobpos(pr, "D2M336"), fake_bc$pheno[,1],
                      contrasts=contrasts)
    expect_equivalent(co1[,"est"], out_fit1a$coef)
    expect_equivalent(co1[,"SE"], out_fit1a$SE)

    out_fit1b <- fit1(pull_genoprobpos(pr, "D5M394"), fake_bc$pheno[,1],
                      contrasts=contrasts)
    expect_equivalent(co2[,"est"], out_fit1b$coef)
    expect_equivalent(co2[,"SE"], out_fit1b$SE)

    # compare to scan1coef()
    ests_c2 <- scan1coef(pr[,"2"], fake_bc$pheno[,1], contrasts=contrasts, se=TRUE)
    expect_equal(ests_c2["D2M336",], out_fit1a$coef)
    expect_equal(attr(ests_c2, "SE")["D2M336",], out_fit1a$SE)

    ests_c5 <- scan1coef(pr[,"5"], fake_bc$pheno[,1], contrasts=contrasts, se=TRUE)
    expect_equal(ests_c5["D5M394",], out_fit1b$coef)
    expect_equal(attr(ests_c5, "SE")["D5M394",], out_fit1b$SE)

})

test_that("fit1 handles contrasts properly in an intercross", {

    library(qtl)
    data(fake.f2)
    fake.f2 <- subset(fake.f2, chr=c(1,13))

    fake_f2 <- convert2cross2(fake.f2)

    fake.f2 <- calc.genoprob(fake.f2, error.prob=0.002, map.function="c-f")
    pr <- calc_genoprob(fake_f2, error_prob=0.002, map_function="c-f")

    # fit single-QTL model at two positions
    qtl <- makeqtl(fake.f2, c(1,13), pos=c(37.1, 24.0), what="prob")
    out.fq1 <- fitqtl(fake.f2, qtl=qtl, formula=y~Q1, get.ests=TRUE, method="hk")
    co1 <- summary(out.fq1)$ests
    out.fq2 <- fitqtl(fake.f2, qtl=qtl, formula=y~Q2, get.ests=TRUE, method="hk")
    co2 <- summary(out.fq2)$ests

    # R/qtl1 and R/qtl2 give different answers for the intercept,
    #    so we'll just look at the additive and dominance effects
    contrasts <- cbind(mu=c(1,1,1), a=c(-1,0,1), d=c(0,1,0))
    out_fit1a <- fit1(pull_genoprobpos(pr, "D1M437"), fake_f2$pheno[,1],
                      contrasts=contrasts)
    expect_equivalent(co1[-1,"est"], out_fit1a$coef[-1])
    expect_equivalent(co1[-1,"SE"], out_fit1a$SE[-1])

    out_fit1b <- fit1(pull_genoprobpos(pr, "D13M254"), fake_f2$pheno[,1],
                      contrasts=contrasts)
    expect_equivalent(co2[-1,"est"], out_fit1b$coef[-1])
    expect_equivalent(co2[-1,"SE"], out_fit1b$SE[-1])

    # compare to scan1coef()
    ests_c1 <- scan1coef(pr[,"1"], fake_f2$pheno[,1], contrasts=contrasts, se=TRUE)
    expect_equal(ests_c1["D1M437",], out_fit1a$coef)
    expect_equal(attr(ests_c1, "SE")["D1M437",], out_fit1a$SE)

    ests_c13 <- scan1coef(pr[,"13"], fake_f2$pheno[,1], contrasts=contrasts, se=TRUE)
    expect_equal(ests_c13["D13M254",], out_fit1b$coef)
    expect_equal(attr(ests_c13, "SE")["D13M254",], out_fit1b$SE)

})

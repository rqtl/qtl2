context("fit1 for binary traits")

test_that("fit1 for binary traits works in intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(18:19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    pheno <- setNames(as.numeric(iron$pheno[,1] > median(iron$pheno[,1])),
                      ind_ids_pheno(iron))
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar, model="binary")

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, model="binary", zerosum=FALSE) })

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(3, 4, 12)
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, model="binary", zerosum=FALSE) })

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
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar, model="binary")

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, se=TRUE, model="binary", zerosum=FALSE) })

    # fit1, missing data
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, se=TRUE, model="binary", zerosum=FALSE) })

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
    glm0 <- glm(pheno ~ covar, family=binomial(link=logit))
    X <- cbind(probs[[1]][,,pmar[1]], covar)
    colnames(X) <- c("SS", "SB", "BB", "ac1")
    glm1 <- glm(pheno ~ -1 + X, family=binomial(link=logit))
    glm_lod <- (glm1$deviance - glm0$deviance)/(-2*log(10))
    p1 <- glm1$fitted
    p0 <- glm0$fitted
    y <- pheno[!is.na(pheno)]
    glm_ind_lod <- (y * log10(p1) + (1-y)*log10(1-p1)) -
        (y * log10(p0) + (1-y)*log10(1-p0))

    expect_equal(out_fit1[[1]]$lod, glm_lod)
    expect_equal(out_fit1[[1]]$ind_lod, glm_ind_lod)

    expect_equal(out_fit1[[1]]$coef, stats::setNames(glm1$coef, c("SS", "SB", "BB", "ac1")))
    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(glm1)$coef[,2], c("SS", "SB", "BB", "ac1")), tol=1e-6)

    # direct calculations, chr X
    glm0 <- glm(pheno ~ Xcovar, family=binomial(link=logit))
    X <- probs[[3]][,,pmar[3]]
    colnames(X) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    glm1 <- glm(pheno ~ -1 + X, family=binomial(link=logit))
    glm_lod <- (glm1$deviance - glm0$deviance)/(-2*log(10))
    p1 <- glm1$fitted
    p0 <- glm0$fitted
    y <- pheno[!is.na(pheno)]
    glm_ind_lod <- (y * log10(p1) + (1-y)*log10(1-p1)) -
        (y * log10(p0) + (1-y)*log10(1-p0))

    expect_equal(out_fit1[[3]]$lod, glm_lod)
    expect_equal(out_fit1[[3]]$ind_lod, glm_ind_lod)
    expect_equal(out_fit1[[3]]$coef, stats::setNames(glm1$coef, c("SS", "SB", "BS", "BB", "SY", "BY")))
    expect_equal(out_fit1[[3]]$SE, stats::setNames(summary(glm1)$coef[,2], c("SS", "SB", "BS", "BB", "SY", "BY")), tol=1e-6)

})


test_that("fit1 by H-K works in riself", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[,4:5]
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)

    pheno <- setNames(as.numeric(grav2$pheno[,219] > 101.25), ind_ids_pheno(grav2))

    # calculate LOD scores
    out <- scan1(probs, pheno, model="binary")

    # estimate coefficients
    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno, model="binary", zerosum=FALSE))

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(1, 172)
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno, model="binary", zerosum=FALSE))

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
    out <- scan1(probs, pheno, model="binary")

    # estimate coefficients
    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno, se=TRUE, model="binary", zerosum=FALSE))

    # fit1, missing data
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno, se=TRUE, model="binary", zerosum=FALSE))

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
    glm0 <- glm(pheno ~ 1, family=binomial(link=logit))
    X <- probs[[1]][,,pmar[1]]
    colnames(X) <- c("LL", "CC")
    glm1 <- glm(pheno ~ -1 + X, family=binomial(link=logit))
    glm_lod <- (glm1$deviance - glm0$deviance)/(-2*log(10))
    p1 <- glm1$fitted
    p0 <- glm0$fitted
    y <- pheno[!is.na(pheno)]
    glm_ind_lod <- (y * log10(p1) + (1-y)*log10(1-p1)) -
        (y * log10(p0) + (1-y)*log10(1-p0))

    expect_equal(out_fit1[[1]]$lod, glm_lod)
    expect_equal(out_fit1[[1]]$ind_lod, glm_ind_lod)
    expect_equal(out_fit1[[1]]$coef, stats::setNames(glm1$coef, c("LL", "CC")))
    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(glm1)$coef[,2], c("LL", "CC")), tol=1e-6)

})

test_that("fit1 works for binary traits with weights", {

    set.seed(17262911)

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(2,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.7)),
                    ind_ids(iron))
    phe[c(108,142,268)] <- NA
    weights <- setNames(sample(1:10, n_ind(iron), replace=TRUE), names(phe))

    npos <- dim(probs)[3,]
    pos <- sapply(npos, sample, 1)
    pr <- list(probs[[1]][,,pos[1]],
               probs[[2]][,,pos[2]])

    out_fit1_1 <- fit1(pr[[1]], phe, model="binary", se=TRUE, weights=weights, zerosum=FALSE)
    out_fit1_2 <- fit1(pr[[2]], phe, model="binary", se=TRUE, weights=weights, zerosum=FALSE)

    # coefficients and SEs
    co2 <- scan1coef(probs[,"2"], phe, model="binary", se=TRUE, weights=weights, zerosum=FALSE)
    expect_equal(out_fit1_1$coef, co2[pos[1],])
    expect_equal(out_fit1_1$SE, attr(co2, "SE")[pos[1],])
    coX <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE, weights=weights, zerosum=FALSE)
    expect_equal(out_fit1_2$coef, coX[pos[2],])
    expect_equal(out_fit1_2$SE, attr(coX, "SE")[pos[2],])

    # lod
    out <- scan1(probs, phe, model="binary", weights=weights)
    expect_equal(out_fit1_1$lod, out[pos[1]])
    expect_equal(out_fit1_2$lod, out[npos[1] + pos[2]])

    # compare to glm
    out_glm_1 <- glm(phe ~ -1 + pr[[1]], family=binomial(link=logit), weights=weights,
                     control=list(epsilon=1e-11))
    out_glm_2 <- glm(phe ~ -1 + pr[[2]], family=binomial(link=logit), weights=weights,
                     control=list(epsilon=1e-11))
    out_glm_0 <- glm(phe ~ 1, family=binomial(link=logit), weights=weights)
    sum_glm_1 <- summary(out_glm_1)
    sum_glm_2 <- summary(out_glm_2)

    # function to calc lod
    lod_glm <- function(out_glm_alt, out_glm_null)
        (out_glm_null$deviance - out_glm_alt$deviance)/(2*log(10))

    # lod scores
    expect_equivalent(out_fit1_1$lod, lod_glm(out_glm_1, out_glm_0))
    expect_equivalent(out_fit1_2$lod, lod_glm(out_glm_2, out_glm_0))

    # coefficients and SEs
    expect_equivalent(out_fit1_1$coef, out_glm_1$coef)
    expect_equivalent(out_fit1_2$coef, out_glm_2$coef)
    expect_equivalent(out_fit1_1$SE, sum_glm_1$coef[,2])
    expect_equivalent(out_fit1_2$SE, sum_glm_2$coef[,2])

    # add a covariate
    X <- setNames(rnorm(n_ind(iron)), names(phe))

    out_fit1_1 <- fit1(pr[[1]], phe, model="binary", se=TRUE, weights=weights, addcovar=X, zerosum=FALSE)
    out_fit1_2 <- fit1(pr[[2]], phe, model="binary", se=TRUE, weights=weights, addcovar=X, zerosum=FALSE)

    # coefficients and SEs
    co2 <- scan1coef(probs[,"2"], phe, model="binary", se=TRUE, weights=weights, addcovar=X, zerosum=FALSE)
    expect_equal(out_fit1_1$coef, co2[pos[1],])
    expect_equal(out_fit1_1$SE, attr(co2, "SE")[pos[1],])
    coX <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE, weights=weights, addcovar=X, zerosum=FALSE)
    expect_equal(out_fit1_2$coef, coX[pos[2],])
    expect_equal(out_fit1_2$SE, attr(coX, "SE")[pos[2],])

    # lod
    out <- scan1(probs, phe, model="binary", weights=weights, addcovar=X)
    expect_equal(out_fit1_1$lod, out[pos[1]])
    expect_equal(out_fit1_2$lod, out[npos[1] + pos[2]])

    # compare to glm
    out_glm_1 <- glm(phe ~ -1 + pr[[1]] + X, family=binomial(link=logit), weights=weights,
                     control=list(epsilon=1e-11))
    out_glm_2 <- glm(phe ~ -1 + pr[[2]] + X, family=binomial(link=logit), weights=weights,
                     control=list(epsilon=1e-11))
    out_glm_0 <- glm(phe ~ X, family=binomial(link=logit), weights=weights)
    sum_glm_1 <- summary(out_glm_1)
    sum_glm_2 <- summary(out_glm_2)

    # lod scores
    expect_equivalent(out_fit1_1$lod, lod_glm(out_glm_1, out_glm_0))
    expect_equivalent(out_fit1_2$lod, lod_glm(out_glm_2, out_glm_0))

    # coefficients and SEs
    expect_equivalent(out_fit1_1$coef, out_glm_1$coef)
    expect_equivalent(out_fit1_2$coef, out_glm_2$coef)
    expect_equivalent(out_fit1_1$SE, sum_glm_1$coef[,2])
    expect_equivalent(out_fit1_2$SE, sum_glm_2$coef[,2])

    # interactive covariate, autosome only
    out_fit1_1 <- fit1(pr[[1]], phe, model="binary", se=TRUE, weights=weights, addcovar=X, intcovar=X, zerosum=FALSE)

    # coefficients and SEs
    co2 <- scan1coef(probs[,"2"], phe, model="binary", se=TRUE, weights=weights, addcovar=X, intcovar=X, zerosum=FALSE)
    expect_equal(out_fit1_1$coef, co2[pos[1],])
    expect_equal(out_fit1_1$SE, attr(co2, "SE")[pos[1],])

    # lod
    out <- scan1(probs, phe, model="binary", weights=weights, addcovar=X, intcovar=X)
    expect_equal(out_fit1_1$lod, out[pos[1]])

    # compare to glm
    out_glm_1 <- glm(phe ~ -1 + pr[[1]] + X + I(pr[[1]][,2]*X) + I(pr[[1]][,3]*X),
                     family=binomial(link=logit), weights=weights,
                     control=list(epsilon=1e-11))
    out_glm_0 <- glm(phe ~ X, family=binomial(link=logit), weights=weights)
    sum_glm_1 <- summary(out_glm_1)

    # lod scores
    expect_equivalent(out_fit1_1$lod, lod_glm(out_glm_1, out_glm_0))

    # coefficients and SEs
    expect_equivalent(out_fit1_1$coef, out_glm_1$coef)
    expect_equivalent(out_fit1_1$SE, sum_glm_1$coef[,2])

})

test_that("fit one for binary traits handles NA case", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c("2", "X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.95)),
                    ind_ids(iron))

    suppressWarnings(out <- scan1(probs, phe, model="binary"))

    mar <- c("D2Mit304", "DXMit186")
    p1 <- pull_genoprobpos(probs, mar[1])
    p2 <- pull_genoprobpos(probs, mar[2])

    suppressWarnings(ll1a <- glm(phe ~ -1 + p1, family=binomial(link=logit))$deviance)
    suppressWarnings(ll1b <- glm(phe ~ -1 + p2, family=binomial(link=logit))$deviance)
    ll0 <- glm(phe ~ 1, family=binomial(link=logit))$deviance

    expect_equal(out[mar[1],1], (ll0-ll1a)/2/log(10), tol=1e-5)
    expect_equal(out[mar[2],1], (ll0-ll1b)/2/log(10), tol=1e-6)

    # repeat with weights
    set.seed(20180720)
    w <- setNames( runif(length(phe), 1, 5), names(phe))
    suppressWarnings(outw <- scan1(probs, phe, model="binary", weights=w))

    suppressWarnings(ll1a_w <- glm(phe ~ -1 + p1, family=binomial(link=logit), weights=w)$deviance)
    suppressWarnings(ll1b_w <- glm(phe ~ -1 + p2, family=binomial(link=logit), weights=w)$deviance)
    suppressWarnings(ll0_w <- glm(phe ~ 1, family=binomial(link=logit), weights=w)$deviance)

    expect_equal(outw[mar[1],1], (ll0_w-ll1a_w)/2/log(10), tol=1e-4)
    expect_equal(outw[mar[2],1], (ll0_w-ll1b_w)/2/log(10), tol=1e-6)

})

test_that("fit one for binary traits handles NA case with DO data", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex, error_prob=0.002, cores=2)
    apr <- genoprob_to_alleleprob(probs)

    phe <- setNames((DOex$pheno[,1] > quantile(DOex$pheno[,1], 0.95, na.rm=TRUE))*1, rownames(DOex$pheno))

    suppressWarnings(out <- scan1(apr, phe, model="binary"))

    mar <- c("backupUNC021331957", "UNC020459944")
    p1 <- pull_genoprobpos(apr, mar[1])
    p2 <- pull_genoprobpos(apr, mar[2])

    suppressWarnings(ll1a <- glm(phe ~ -1 + p1, family=binomial(link=logit))$deviance)
    suppressWarnings(ll1b <- glm(phe ~ -1 + p2, family=binomial(link=logit))$deviance)
    ll0 <- glm(phe ~ 1, family=binomial(link=logit))$deviance

    expect_equal(out[mar[1],1], (ll0-ll1a)/2/log(10))
    expect_equal(out[mar[2],1], (ll0-ll1b)/2/log(10))


    # repeat with weights
    set.seed(20180720)
    w <- setNames( runif(length(phe), 1, 5), names(phe))
    suppressWarnings(outw <- scan1(apr, phe, model="binary", weights=w))

    suppressWarnings(ll1a_w <- glm(phe ~ -1 + p1, family=binomial(link=logit), weights=w)$deviance)
    suppressWarnings(ll1b_w <- glm(phe ~ -1 + p2, family=binomial(link=logit), weights=w)$deviance)
    suppressWarnings(ll0_w <- glm(phe ~ 1, family=binomial(link=logit), weights=w)$deviance)

    expect_equal(outw[mar[1],1], (ll0_w-ll1a_w)/2/log(10))
    expect_equal(outw[mar[2],1], (ll0_w-ll1b_w)/2/log(10), tol=1e-6)

})

context("fit1 for binary traits")

test_that("fit1 for binary traits works in intercross", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
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
### FIX_ME: scan1coef for binary traits not yet implemented
#    coef <- lapply(seq_len(length(probs)), function(i) {
#        if(i==3) cov <- NULL
#        else cov <- covar
#        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, model="binary") })

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(3, 4, 12)
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov, model="binary") })

    pos <- cumsum(c(0, npos[-3])) + pmar
    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:3) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:3)
#        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # repeat the whole thing with a couple of missing phenotypes
    pheno[c(187, 244)] <- NA

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients; no covariates for X chromosome
### FIX_ME: scan1coef for binary traits not yet implemented
#    coef <- lapply(seq_len(length(probs)), function(i) {
#        if(i==3) cov <- NULL
#        else cov <- covar
#        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov, se=TRUE) })

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
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:3)
#        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # check SEs
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:3)
#        expect_equal(out_fit1[[i]]$SE, attr(coef[[i]], "SE")[pmar[i],])

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

    expect_equal(out_fit1[[1]]$lod, glm_lod, tol=1e-3) ## FIX_ME this just looks wrong
    expect_equal(out_fit1[[1]]$ind_lod, glm_ind_lod, tol=1e-3) ## FIX_ME this just looks wrong

### FIX_ME: scan1coef for binary traits not yet implemented
#    expect_equal(out_fit1[[1]]$coef, stats::setNames(lm1$coef, c("SS", "SB", "BB", "ac1")))
#    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(lm1)$coef[,2], c("SS", "SB", "BB", "ac1")))

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

    expect_equal(out_fit1[[3]]$lod, glm_lod, tol=0.015) ## FIX_ME this just looks wrong
    expect_equal(out_fit1[[3]]$ind_lod, glm_ind_lod, tol=0.0086) ## FIX_ME this just looks wrong
### FIX_ME: scan1coef for binary traits not yet implemented
#    expect_equal(out_fit1[[3]]$coef, stats::setNames(lm1$coef, c("SS", "SB", "BS", "BB", "SY", "BY")))
#    expect_equal(out_fit1[[3]]$SE, stats::setNames(summary(lm1)$coef[,2], c("SS", "SB", "BS", "BB", "SY", "BY")))

})


test_that("fit1 by H-K works in riself", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    grav2 <- grav2[,4:5]
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)

    pheno <- setNames(as.numeric(grav2$pheno[,219] > 101.25), ind_ids_pheno(grav2))

    # calculate LOD scores
    out <- scan1(probs, pheno, model="binary")

    # estimate coefficients
### FIX_ME: scan1coef for binary traits not yet implemented
#    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno))

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(1, 172)
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno, model="binary"))

    pos <- c(0,npos[1]) + pmar
    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:2) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:2)
#        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # repeat the whole thing with a couple of missing phenotypes
    pheno[c(24, 106)] <- NA

    # calculate LOD scores
    out <- scan1(probs, pheno, model="binary")

    # estimate coefficients
### FIX_ME: scan1coef for binary traits not yet implemented
#    coef <- lapply(seq_len(length(probs)), function(i) scan1coef(subset(probs, chr=names(probs)[i]), pheno, se=TRUE))

    # fit1, no missing data
    out_fit1 <- lapply(seq(along=pmar), function(i) fit1(probs[[i]][,,pmar[i]], pheno, se=TRUE, model="binary"))

    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:2) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:2)
#        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])


    # check SEs
### FIX_ME: scan1coef for binary traits not yet implemented
#    for(i in 1:2)
#        expect_equal(out_fit1[[i]]$SE, attr(coef[[i]], "SE")[pmar[i],])

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
### FIX_ME: scan1coef for binary traits not yet implemented
#    expect_equal(out_fit1[[1]]$coef, stats::setNames(lm1$coef, c("LL", "CC")))
#    expect_equal(out_fit1[[1]]$SE, stats::setNames(summary(lm1)$coef[,2], c("LL", "CC")))

})

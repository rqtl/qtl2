context("scan1coef for binary traits")

test_that("scan1coef for binary traits works with intercross", {

    set.seed(17262911)

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(2,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.7)),
                    ind_ids(iron))

    # no covariates, autosome
    co <- scan1coef(probs[,"2"], phe, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # no covariates, X chromosome
    co <- scan1coef(probs[,"X"], phe, model="binary")
    coSE <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    sex <- setNames(as.numeric(iron$is_female), names(iron$is_female))

    # add covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a + sex, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    covar <- setNames(rnorm(n_ind(iron)), names(sex))

    # add covariate, X chromosome
    co <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary")
    coSE <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a + covar, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # int covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) {
        X <- cbind(a, sex, a[,-1]*sex)
        glm(phe ~ -1 + X, family=binomial(link=logit),
            control=list(epsilon=1e-12)) })
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1", "SB:ic1", "BB:ic1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

})


test_that("scan1coef for binary traits works some missing phenotypes", {

    set.seed(17262911)

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(2,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.7)),
                    ind_ids(iron))
    phe[c(41,153)] <- NA

    # no covariates, autosome
    co <- scan1coef(probs[,"2"], phe, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # no covariates, X chromosome
    co <- scan1coef(probs[,"X"], phe, model="binary")
    coSE <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    sex <- setNames(as.numeric(iron$is_female), names(iron$is_female))

    # add covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a + sex, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    covar <- setNames(rnorm(n_ind(iron)), names(sex))

    # add covariate, X chromosome
    co <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary")
    coSE <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a + covar, family=binomial(link=logit),
                                                    control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # int covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary")
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) {
        X <- cbind(a, sex, a[,-1]*sex)
        glm(phe ~ -1 + X, family=binomial(link=logit),
            control=list(epsilon=1e-12)) })
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1", "SB:ic1", "BB:ic1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

})

test_that("scan1coef for binary traits works with weights", {

    set.seed(17262911)

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(2,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.7)),
                    ind_ids(iron))
    weights <- setNames(sample(1:10, n_ind(iron), replace=TRUE), names(phe))

    # no covariates, autosome
    co <- scan1coef(probs[,"2"], phe, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, model="binary", , weights=weights, se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # no covariates, X chromosome
    co <- scan1coef(probs[,"X"], phe, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    sex <- setNames(as.numeric(iron$is_female), names(iron$is_female))

    # add covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a + sex, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    covar <- setNames(rnorm(n_ind(iron)), names(sex))

    # add covariate, X chromosome
    co <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a + covar, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # int covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) {
        X <- cbind(a, sex, a[,-1]*sex)
        glm(phe ~ -1 + X, family=binomial(link=logit),
            control=list(epsilon=1e-12), weights=weights) })
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1", "SB:ic1", "BB:ic1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

})

test_that("scan1coef for binary traits works with weights and missing phenotypes", {

    set.seed(17262911)

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(2,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, err=0.002)

    phe <- iron$pheno[,1]
    phe <- setNames(as.numeric(phe > quantile(phe, 0.7)),
                    ind_ids(iron))
    phe[c(108,142,268)] <- NA
    weights <- setNames(sample(1:10, n_ind(iron), replace=TRUE), names(phe))

    # no covariates, autosome
    co <- scan1coef(probs[,"2"], phe, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, model="binary", , weights=weights, se=TRUE)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # no covariates, X chromosome
    co <- scan1coef(probs[,"X"], phe, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"X"], phe, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    sex <- setNames(as.numeric(iron$is_female), names(iron$is_female))

    # add covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) glm(phe ~ -1 + a + sex, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    covar <- setNames(rnorm(n_ind(iron)), names(sex))

    # add covariate, X chromosome
    co <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"X"], phe, addcovar=covar, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[2]], 3, function(a) glm(phe ~ -1 + a + covar, family=binomial(link=logit),
                                                    weights=weights, control=list(epsilon=1e-12)))
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BS", "BB", "SY", "BY", "ac1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

    # int covariate, autosome
    co <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", weights=weights)
    coSE <- scan1coef(probs[,"2"], phe, addcovar=sex, intcovar=sex, model="binary", se=TRUE, weights=weights)
    expect_equivalent(co, coSE)

    out_glm <- apply(probs[[1]], 3, function(a) {
        X <- cbind(a, sex, a[,-1]*sex)
        glm(phe ~ -1 + X, family=binomial(link=logit),
            control=list(epsilon=1e-12), weights=weights) })
    glm_coef <- t(sapply(out_glm, function(a) a$coef))
    glm_se <- t(sapply(out_glm, function(a) summary(a)$coef[,2]))
    colnames(glm_se) <- colnames(glm_coef) <- c("SS", "SB", "BB", "ac1", "SB:ic1", "BB:ic1")
    expect_equivalent(unclass(co), glm_coef)
    expect_equal(attr(coSE, "SE"), glm_se, tol=1e-6)

})

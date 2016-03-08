context("effect scan by scan1coef")

convert_probs2qtl2 <-
    function(cross)
{
    list(probs=lapply(cross$geno, function(a) {
        pr <- aperm(a$prob, c(1,3,2))
        rownames(pr) <- paste(1:nrow(pr))
        pr }),
         map=lapply(cross$geno, function(a) attr(a$prob, "map")),
         indID=paste(1:qtl::nind(cross)),
         chrID=names(cross$geno))

}

# calc estimates via lm(), just one chromosome
eff_via_lm <-
    function(probs, pheno, addcovar=NULL, intcovar=NULL, weights=NULL,
             se=TRUE)
{
    npos <- dim(probs)[2]
    nind <- length(pheno)
    ngen <- dim(probs)[3]

    nadd <- ifelse(is.null(addcovar), 0, ncol(addcovar))
    nint <- ifelse(is.null(intcovar), 0, ncol(intcovar))
    ncoef <- ngen + nadd + (ngen-1)*nint

    result <- matrix(nrow=npos, ncol=ncoef)
    SEs <- matrix(nrow=npos, ncol=ncoef)

    if(is.null(weights)) wts <- rep(1, nind)
    else wts <- weights
    if(!is.null(intcovar)) intcovar <- as.matrix(intcovar)

    for(i in 1:npos) {
        X <- probs[,i,]
        if(!is.null(addcovar)) X <- cbind(X, addcovar)
        if(!is.null(intcovar)) {
            for(j in 1:ncol(intcovar))
                X <- cbind(X, probs[,i,-1]*intcovar[,j])
        }

        lm.out <- lm(pheno ~ -1 + X, weights=weights)
        result[i,] <- lm.out$coef
        if(se) # need to deal with NAs
            SEs[i,!is.na(lm.out$coef)] <- summary(lm.out)$coef[,2]
    }

    if(se) attr(result, "SE") <- SEs

    class(result) <- c("scan1coef", "scan1", "matrix")
    result
}

test_that("scan1coef for backcross", {

    set.seed(9308594)

    library(qtl)
    data(hyper)
    hyper <- hyper[4,] # chr 4 only
    hyper <- calc.genoprob(hyper, step=5)
    prob <- hyper$geno[[1]]$prob
    phe <- hyper$pheno[,1]
    weights <- runif(nind(hyper), 0, 5)
    covar <- cbind(sex=sample(0:1, nind(hyper), replace=TRUE))
    rownames(prob) <- names(phe) <- names(weights) <- rownames(covar) <- paste(1:nind(hyper))

    # probs for qtl2scan code
    prob2 <- convert_probs2qtl2(hyper)

    # no covariates
    expected <- eff_via_lm(prob, phe)
    coef <- scan1coef(prob2, phe)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA")))
    coef <- scan1coef(prob2, phe, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA")))

    # no covariates, weighted
    expected <- eff_via_lm(prob, phe, weights=weights)
    coef <- scan1coef(prob2, phe, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA")))
    coef <- scan1coef(prob2, phe, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA")))

    # one add've covariate
    expected <- eff_via_lm(prob, phe, covar)
    coef <- scan1coef(prob2, phe, covar)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex")))
    coef <- scan1coef(prob2, phe, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA", "sex")))

    # one add've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex")))
    coef <- scan1coef(prob2, phe, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA", "sex")))

    # one int've covariate
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(colnames(prob), c("BB", "BA", "sex", "BA:sex")))

    # two int've covariates
    covar <- cbind(covar, another=rnorm(nind(hyper)))
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))

    # two int've covariates, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(colnames(prob), c("BB", "BA", "sex", "another", "BA:sex", "BA:another")))

})

test_that("scan1coef for backcross, with contrasts", {

    set.seed(9308594)

    library(qtl)
    data(hyper)
    hyper <- hyper[4,] # chr 4 only
    hyper <- calc.genoprob(hyper, step=5)
    prob <- hyper$geno[[1]]$prob
    phe <- hyper$pheno[,1]
    weights <- runif(nind(hyper), 0, 5)
    covar <- cbind(sex=sample(0:1, nind(hyper), replace=TRUE))
    rownames(prob) <- names(phe) <- names(weights) <- rownames(covar) <- paste(1:nind(hyper))

    # probs for qtl2scan code
    prob2 <- convert_probs2qtl2(hyper)

    # use contrasts
    contrasts <- cbind(mu=c(1,1), a=c(-0.5, 0.5))

    # change to contrasts for use with lm()
    p <- prob
    dim(p) <- c(prod(dim(prob)[1:2]), dim(prob)[3])
    p <- p %*% contrasts
    dim(p) <- dim(prob)
    prob <- p
    dimnames(prob) <- list(rownames(prob2), dimnames(prob2)[[3]], colnames(contrasts))

    posnames <- dimnames(prob2$probs[[1]])[[3]]

    # no covariates
    expected <- eff_via_lm(prob, phe)
    coef <- scan1coef(prob2, phe, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a")))
    coef <- scan1coef(prob2, phe, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a")))

    # no covariates, weighted
    expected <- eff_via_lm(prob, phe, weights=weights)
    coef <- scan1coef(prob2, phe, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a")))
    coef <- scan1coef(prob2, phe, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a")))

    # one add've covariate
    expected <- eff_via_lm(prob, phe, covar)
    coef <- scan1coef(prob2, phe, covar, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex")))
    coef <- scan1coef(prob2, phe, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "sex")))

    # one add've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex")))
    coef <- scan1coef(prob2, phe, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "sex")))

    # one int've covariate
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "a:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "a:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "sex", "a:sex")))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "a:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "a:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "sex", "a:sex")))

    # two int've covariates
    covar <- cbind(covar, another=rnorm(nind(hyper)))
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))

    # two int've covariates, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(posnames, c("mu", "a", "sex", "another", "a:sex", "a:another")))

})



test_that("scan1coef for intercross", {

    set.seed(9308594)

    library(qtl2scan)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    prob2 <- calc_genoprob(iron[,7], step=1, error_prob=0.002) # just look at chr 7
    phe <- iron$pheno[,1] # liver phenotype
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    covar <- cbind(sex=covar)
    weights <- runif(n_ind(iron), 0, 5)
    names(weights) <- names(phe)

    # different organization of probs for qtl2scan and lm() code
    prob <- aperm(prob2$probs[[1]], c(1,3,2)) # rearrange as expected for lm()

    # no covariates
    expected <- eff_via_lm(prob, phe)
    coef <- scan1coef(prob2, phe)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    posnames <- dimnames(prob2$probs[[1]])[[3]]

    # no covariates, weighted
    expected <- eff_via_lm(prob, phe, weights=weights)
    coef <- scan1coef(prob2, phe, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB")))
    coef <- scan1coef(prob2, phe, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("SS", "SB", "BB")))

    # one add've covariate
    expected <- eff_via_lm(prob, phe, covar)
    coef <- scan1coef(prob2, phe, covar)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex")))
    coef <- scan1coef(prob2, phe, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("SS", "SB", "BB", "sex")))

    # one add've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex")))
    coef <- scan1coef(prob2, phe, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("SS", "SB", "BB", "sex")))

    # one int've covariate
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("SS", "SB", "BB", "sex", "SB:sex", "BB:sex")))

    # two int've covariates
    covar <- cbind(covar, another=rnorm(n_ind(iron)))
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equal(dimnames(coef),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))

    # two int've covariates, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")),
                 list(posnames, c("SS", "SB", "BB", "sex", "another", "SB:sex", "BB:sex",
                                        "SB:another", "BB:another")))

})

test_that("scan1coef for intercross, with contrasts", {

    set.seed(9308594)

    library(qtl2scan)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    prob2 <- calc_genoprob(iron[,7], step=1, error_prob=0.002) # just look at chr 7
    phe <- iron$pheno[,1] # liver phenotype
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    covar <- cbind(sex=covar)
    weights <- runif(n_ind(iron), 0, 5)
    names(weights) <- names(phe)

    # different organization of probs for qtl2scan and lm() code
    prob <- aperm(prob2$probs[[1]], c(1,3,2)) # rearrange as expected for lm()

    # use contrasts
    contrasts <- cbind(mu=c(1,1,1), a=c(-0.5, 0, 0.5), d=c(-0.5, 1, -0.5))

    # change to contrasts for use with lm()
    dn <- dimnames(prob)
    p <- prob
    dim(p) <- c(prod(dim(prob)[1:2]), dim(prob)[3])
    p <- p %*% contrasts
    dim(p) <- dim(prob)
    prob <- p
    dimnames(prob) <- list(dn[[1]], dn[[2]], colnames(contrasts))

    posnames <- dimnames(prob2$probs[[1]])[[3]]

    # no covariates
    expected <- eff_via_lm(prob, phe)
    coef <- scan1coef(prob2, phe, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d")))
    coef <- scan1coef(prob2, phe, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d")))

    # no covariates, weighted
    expected <- eff_via_lm(prob, phe, weights=weights)
    coef <- scan1coef(prob2, phe, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d")))
    coef <- scan1coef(prob2, phe, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d")))

    # one add've covariate
    expected <- eff_via_lm(prob, phe, covar)
    coef <- scan1coef(prob2, phe, covar, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex")))
    coef <- scan1coef(prob2, phe, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex")))

    # one add've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex")))
    coef <- scan1coef(prob2, phe, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex")))

    # one int've covariate
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex", "a:sex", "d:sex")))

    # two int've covariates
    covar <- cbind(covar, another=rnorm(n_ind(iron)))
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "another",
                                                              "a:sex", "d:sex", "a:another", "d:another")))
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "another",
                                                              "a:sex", "d:sex", "a:another", "d:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex", "another",
                                                              "a:sex", "d:sex", "a:another", "d:another")))

    # two int've covariates, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames,
                                      c("mu", "a", "d", "sex", "another", "a:sex", "d:sex",
                                        "a:another", "d:another")))
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, contrasts=contrasts, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equal(dimnames(coef), list(posnames, c("mu", "a", "d", "sex", "another",
                                                              "a:sex", "d:sex", "a:another", "d:another")))
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))
    expect_equal(dimnames(attr(coef, "SE")), list(posnames, c("mu", "a", "d", "sex", "another",
                                                              "a:sex", "d:sex", "a:another", "d:another")))

})

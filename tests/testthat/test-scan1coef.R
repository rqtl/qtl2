context("effect scan by scan1coef")

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
                X <- cbind(X, probs[,i,-1,drop=FALSE]*intcovar[,j])
        }

        lm.out <- lm(pheno ~ -1 + X, weights=weights)
        result[i,] <- lm.out$coef
        if(se) # need to deal with NAs
            SEs[i,!is.na(lm.out$coef)] <- summary(lm.out)$coef[,2]
    }

    if(se) attr(result, "SE") <- SEs

    result
}

test_that("scan1coef for backcross", {

    set.seed(9308594)

    library(qtl)
    data(hyper)
    hyper <- hyper[4,] # chr 4 only
    prob <- calc.genoprob(hyper, step=5)$geno[[1]]$prob
    phe <- hyper$pheno[,1]
    weights <- runif(nind(hyper), 0, 5)
    covar <- cbind(sex=sample(0:1, nind(hyper), replace=TRUE))
    rownames(prob) <- names(phe) <- names(weights) <- rownames(covar) <- paste(1:nind(hyper))

    # no covariates
    expected <- eff_via_lm(prob, phe)

    prob2 <- aperm(prob, c(1,3,2)) # rearrange as in R/qtl2
    coef <- scan1coef(prob2, phe)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # no covariates, weighted
    expected <- eff_via_lm(prob, phe, weights=weights)
    coef <- scan1coef(prob2, phe, weights=weights)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # one add've covariate
    expected <- eff_via_lm(prob, phe, covar)
    coef <- scan1coef(prob2, phe, covar)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # one add've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, weights=weights)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # one int've covariate
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # two int've covariate
    covar <- cbind(covar, rnorm(nind(hyper)))
    expected <- eff_via_lm(prob, phe, covar, covar)
    coef <- scan1coef(prob2, phe, covar, covar)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

    # one int've covariate, weighted
    expected <- eff_via_lm(prob, phe, covar, covar, weights=weights)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights)
    expect_equivalent(coef, expected)
    coef <- scan1coef(prob2, phe, covar, covar, weights=weights, se=TRUE)
    expect_equivalent(coef, expected)
    expect_equivalent(attr(coef, "SE"), attr(expected, "SE"))

})

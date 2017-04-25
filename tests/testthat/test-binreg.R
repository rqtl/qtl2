context("logistic regression");

test_that("binreg_eigen functions work", {

    set.seed(94651642)

    n <- 200
    X <- cbind(1, sample(0:1, n, replace=TRUE), rnorm(n))
    nu <- X %*% c(0, 0.5, 0.3)
    p <- exp(nu)/(1+exp(nu))
    y <- rbinom(n, 1, p)

    out <- glm(y ~ -1 + X, family=binomial(link=logit))
    pi <- out$fitted
    expected <- sum(y*log10(pi) + (1-y)*log10(1-pi))

    # just the log likelihood
    expect_equal(calc_ll_binreg_eigenchol(X, y), expected)
    expect_equal(calc_ll_binreg_eigenqr(X, y), expected)
    expect_equal(calc_ll_binreg(X, y), expected)
    out_fit_binreg <- fit_binreg(X, y)
    expect_equal(out_fit_binreg$log10lik, expected);
    expect_equal(expected, out$deviance/(-2*log(10)))


    # coefficients
    coef <- setNames(out$coef, NULL)
    expect_equal(calc_coef_binreg_eigenqr(X, y), coef)
    expect_equal(calc_coef_binreg(X, y), coef)
    out_coefSE_binreg <- calc_coefSE_binreg(X, y)
    expect_equal(out_coefSE_binreg$coef, coef)
    expect_equal(out_fit_binreg$coef, coef)

    # SEs
    SE <- setNames(summary(out)$coef[,2], NULL)
    expect_equal(out_coefSE_binreg$SE, SE, tol=1e-6)
    expect_equal(out_fit_binreg$SE, SE, tol=1e-6)

    # reduced rank matrix
    XX <- cbind(X, X[,1]+X[,2])

    out <- glm(y ~ -1 + XX, family=binomial(link=logit))
    pi <- out$fitted
    expected2 <- sum(y*log10(pi) + (1-y)*log10(1-pi))

    # just the log likelihood
    expect_equal(calc_ll_binreg_eigenchol(XX, y), expected2)
    expect_equal(calc_ll_binreg_eigenqr(XX, y), expected2)
    expect_equal(calc_ll_binreg(XX, y), expected2)
    out_fit_binreg <- fit_binreg(XX, y)
    expect_equal(out_fit_binreg$log10lik, expected2);
    expect_equal(expected, out$deviance/(-2*log(10)))

})


test_that("weighted binreg_eigen functions work", {

    set.seed(94651642)

    n <- 200
    X <- cbind(1, sample(0:1, n, replace=TRUE), rnorm(n))
    nu <- X %*% c(0, 0.5, 0.3)
    p <- exp(nu)/(1+exp(nu))
    wt <- sample(1:10, n, replace=TRUE)
    y <- rbinom(n, wt, p)/wt

    out <- glm(y ~ -1 + X, family=binomial(link=logit), weights=wt)
    pi <- out$fitted
    expected <- sum((y*log10(pi) + (1-y)*log10(1-pi))*wt)

    # just the log likelihood
    expect_equal(calc_ll_binreg_weighted_eigenchol(X, y, wt), expected)
    expect_equal(calc_ll_binreg_weighted_eigenqr(X, y, wt), expected)
    expect_equal(calc_ll_binreg_weighted(X, y, wt), expected)
    out_fit_binreg_weighted <- fit_binreg_weighted(X, y, wt)
    expect_equal(out_fit_binreg_weighted$log10lik, expected);

    # not sure about constant factor in deviance, so calc LOD vs model with just intercept
    out0 <- glm(y ~ 1, family=binomial(link=logit), weights=wt)
    expect_equal(calc_ll_binreg_weighted_eigenqr(X, y, wt) -
                 calc_ll_binreg_weighted_eigenqr(cbind(rep(1, n)), y, wt),
                 (out$deviance - out0$deviance)/(-2*log(10)))

    # coefficients
    coef <- setNames(out$coef, NULL)
    expect_equal(calc_coef_binreg_weighted_eigenqr(X, y, wt), coef)
    expect_equal(calc_coef_binreg_weighted(X, y, wt), coef)
    out_coefSE_binreg_weighted <- calc_coefSE_binreg_weighted(X, y, wt)
    expect_equal(out_coefSE_binreg_weighted$coef, coef)
    expect_equal(out_fit_binreg_weighted$coef, coef)

    # SEs
    SE <- setNames(summary(out)$coef[,2], NULL)
    expect_equal(out_coefSE_binreg_weighted$SE, SE, tol=1e-6)
    expect_equal(out_fit_binreg_weighted$SE, SE, tol=1e-6)

    # reduced rank matrix
    XX <- cbind(X, X[,1]+X[,2])

    out <- glm(y ~ -1 + XX, family=binomial(link=logit), weights=wt)
    pi <- out$fitted
    expected2 <- sum((y*log10(pi) + (1-y)*log10(1-pi))*wt)

    # just the log likelihood
    expect_equal(calc_ll_binreg_weighted_eigenchol(XX, y, wt), expected2)
    expect_equal(calc_ll_binreg_weighted_eigenqr(XX, y, wt), expected2)
    expect_equal(calc_ll_binreg_weighted(XX, y, wt), expected2)
    out_fit_binreg_weighted <- fit_binreg_weighted(XX, y, wt)
    expect_equal(out_fit_binreg_weighted$log10lik, expected2);

    # not sure about constant factor in deviance, so calc LOD vs model with just intercept
    out0 <- glm(y ~ 1, family=binomial(link=logit), weights=wt)
    expect_equal(calc_ll_binreg_weighted_eigenqr(XX, y, wt) -
                 calc_ll_binreg_weighted_eigenqr(cbind(rep(1, n)), y, wt),
                 (out$deviance - out0$deviance)/(-2*log(10)))

})

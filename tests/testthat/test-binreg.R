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

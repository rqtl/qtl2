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

    expect_equal(calc_ll_binreg_eigenchol(X, y), expected)

    # reduced rank matrix
    XX <- cbind(X, X[,1]+X[,2])

    out <- glm(y ~ -1 + XX, family=binomial(link=logit))
    pi <- out$fitted
    expected2 <- sum(y*log10(pi) + (1-y)*log10(1-pi))

    expect_equal(calc_ll_binreg_eigenchol(XX, y), expected2)

})

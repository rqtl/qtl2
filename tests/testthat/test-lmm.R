context("linear mixed models")

set.seed(20151124)
n <- 100
x <- matrix(rnorm(20*n), ncol=20)
k <- as.matrix(-dist(x))
k <- (k-min(k))/(max(k)-min(k))
dimnames(k) <- NULL

X <- cbind(1, sample(0:1, n, replace=TRUE))
d <- chol(0.25*k + diag(rep(0.75, n)))
y <- as.numeric(matrix(rnorm(n), nrow=1) %*% d) + X %*% c(1, 2)


test_that("eigen decomposition works", {

    e <- Rcpp_eigen_decomp(k)

    # same determinants?
    expect_equal(sum(log(e$values)), sum(log(eigen(k)$values)))

    # multiply back
    expect_equal(t(e$vectors) %*% diag(e$values) %*% e$vectors, k)

    # inverse
    expect_equal(t(e$vectors) %*% diag(1/e$values) %*% e$vectors, solve(k))

})

test_that("eigen + rotation works", {

    e <- Rcpp_eigen_rotation(k, y, X)

    # rotate back
    expect_equal(t(e$Kve_t) %*% e$y, y)
    expect_equal(t(e$Kve_t) %*% e$X, X)

    # same determinants?
    expect_equal(sum(log(e$Kva)), sum(log(eigen(k)$values)))

    # multiply back
    expect_equal(t(e$Kve_t) %*% diag(e$Kva) %*% e$Kve_t, k)

    # inverse
    expect_equal(t(e$Kve_t) %*% diag(1/e$Kva) %*% e$Kve_t, solve(k))

})

test_that("fitLMM works", {

    expected_reml <- list(loglik = -217.595158548486,
                          hsq = 0.0367120313458048,
                          sigmasq = 0.79659576062727,
                          beta = c(0.886161208390094, 1.74687862579568))

    expected_ml <- list(loglik = -217.180784548857,
                        hsq = 0,
                        sigmasq = 0.769853921083333,
                        beta = c(0.887111992802203, 1.74121566369042))

    e <- Rcpp_eigen_rotation(k, y, X)
    expect_equal(Rcpp_fitLMM(e$Kva, e$y, e$X), expected_reml)
    expect_equal(Rcpp_fitLMM(e$Kva, e$y, e$X, FALSE), expected_ml)

})

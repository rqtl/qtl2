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

    expected_reml <- structure(list(loglik = -218.467712411475, hsq = 0.0743352779641662,
                                    sigmasq = 0.819543294059229,
                                    beta = c(0.826587328164329, 1.84713422450015)),
                               .Names = c("loglik", "hsq", "sigmasq", "beta"))

    expected_ml <- structure(list(loglik = -217.87551923171, hsq = 0, sigmasq = 0.780625465836864,
                                  beta = c(0.824833750089415, 1.84392110362325)),
                             .Names = c("loglik", "hsq", "sigmasq", "beta"))


    e <- Rcpp_eigen_rotation(k, y, X)
    expect_equal(Rcpp_fitLMM(e$Kva, e$y, e$X), expected_reml)
    expect_equal(Rcpp_fitLMM(e$Kva, e$y, e$X, FALSE), expected_ml)

})

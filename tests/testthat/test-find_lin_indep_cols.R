context("find linearly independent cols")

test_that("find_lin_indep_cols works", {

    set.seed(20151130)
    x <- cbind(1, sample(0:1, 200, repl=TRUE))

    expect_equal(sort(find_lin_indep_cols(x)), c(1,2))
    expect_equal(sort(find_lin_indep_cols(cbind(x, 1))), c(1,2))
    expect_equal(sort(find_lin_indep_cols(cbind(x, x[,1] + 0.5*x[,2]))), c(2,3))

    X <- matrix(rnorm(1000*10), ncol=10)
    expect_equal(sort(find_lin_indep_cols(X)), 1:10)
    expect_equal(sort(find_lin_indep_cols(cbind(rowSums(X), X))), c(1, 3:11))

})

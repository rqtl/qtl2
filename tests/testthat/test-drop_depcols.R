context("drop linearly dependent columns")

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

test_that("drop_depcols works", {

    set.seed(20151130)
    n.ind <- 100
    x <- cbind(1,
               sample(0:1, n.ind, replace=TRUE),
               sample(0:1, n.ind, replace=TRUE))

    expect_equal(drop_depcols(NULL), NULL)
    expect_equal(drop_depcols(x), x)
    for(i in 1:ncol(x))
        expect_equal(drop_depcols(x[,i]), x[,i,drop=FALSE])

    X <- cbind(rowSums(x), x)
    expect_equal(drop_depcols(X), X[,c(1,3,4)])

})

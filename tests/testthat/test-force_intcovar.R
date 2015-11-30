context("force intcovar into addcovar")

test_that("force_intcovar works", {

    set.seed(20151130)
    n.ind <- 100
    X <- cbind(1,
               sample(0:1, n.ind, replace=TRUE),
               sample(0:1, n.ind, replace=TRUE))

    expect_equal(force_intcovar(NULL, NULL), NULL)

    expect_equal(force_intcovar(X, NULL), X)
    expect_equal(force_intcovar(NULL, X), X)

    expect_equal(force_intcovar(X, X), X)

    for(i in 1:ncol(X))
        expect_equal(force_intcovar(X, X[,i]), X)

    expect_equal(force_intcovar(X[,1], X), X)
    expect_equal(force_intcovar(X[,2], X), X[,c(2,1,3)])
    expect_equal(force_intcovar(X[,3], X), X[,c(3,1,2)])

    expect_equal(force_intcovar(X[,-1], X[,1]), X[,c(2,3,1)])
    expect_equal(force_intcovar(X[,-2], X[,2]), X[,c(1,3,2)])
    expect_equal(force_intcovar(X[,-3], X[,3]), X)

})

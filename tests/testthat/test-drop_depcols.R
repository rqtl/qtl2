context("drop linearly dependent columns")

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

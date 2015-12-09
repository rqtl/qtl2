context("force intcovar into addcovar")

test_that("find_matching_cols works", {

    set.seed(20151130)
    x <- cbind(1, sample(0:1, 200, repl=TRUE))

    expect_equal(find_matching_cols(x), c(-1, -1))
    expect_equal(find_matching_cols(cbind(x, 1)), c(-1, -1, 1))
    expect_equal(find_matching_cols(cbind(x, 2, x[,2])), c(-1, -1, -1, 2))
    expect_equal(find_matching_cols(cbind(x, 2, x[,2]+0.01)), c(-1, -1, -1, -1))
    expect_equal(find_matching_cols(cbind(x, 2, 1-x[,2])), c(-1, -1, -1, -1))

    y <- x
    y[5,1] <- NA
    expect_equal(find_matching_cols(cbind(x,y)), c(-1, -1, -1, 2))
    expect_equal(find_matching_cols(cbind(x,y,y[,1])), c(-1, -1, -1, 2, 3))

})


test_that("force_intcovar works", {

    set.seed(20151130)
    n.ind <- 100
    X <- cbind(1,
               sample(0:1, n.ind, replace=TRUE),
               sample(0:1, n.ind, replace=TRUE))
    rownames(X) <- paste(1:n.ind)

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


test_that("force_intcovar works with missing data", {

    set.seed(20151130)
    n.ind <- 100
    X <- cbind(1,
               sample(0:1, n.ind, replace=TRUE),
               sample(0:1, n.ind, replace=TRUE))
    rownames(X) <- paste(1:n.ind)

    # add a few missing values
    X[1,1] <- X[2:3,2] <- X[3:5,3] <- NA

    expect_equal(force_intcovar(X, NULL), X)
    expect_equal(force_intcovar(NULL, X), X)

    expect_equal(force_intcovar(X, X), X[-(1:5),])

    for(i in 1:ncol(X))
        expect_equal(force_intcovar(X, X[,i]), X[-(1:5),])

    expect_equal(force_intcovar(X[,1], X), X[-(1:5),])
    expect_equal(force_intcovar(X[,2], X), X[-(1:5),c(2,1,3)])
    expect_equal(force_intcovar(X[,3], X), X[-(1:5),c(3,1,2)])

    expect_equal(force_intcovar(X[,-1], X[,1]), X[-(1:5),c(2,3,1)])
    expect_equal(force_intcovar(X[,-2], X[,2]), X[-(1:5),c(1,3,2)])
    expect_equal(force_intcovar(X[,-3], X[,3]), X[-(1:5),])

    expect_equal(force_intcovar(X[,1], X[,2]), X[-(1:3),-3])
    expect_equal(force_intcovar(X[,1], X[,3]), X[-c(1,3:5),-2])
    expect_equal(force_intcovar(X[,2], X[,1]), X[-(1:3),c(2,1)])
    expect_equal(force_intcovar(X[,2], X[,3]), X[-(2:5), -1])
    expect_equal(force_intcovar(X[,3], X[,1]), X[-c(1,3:5),c(3,1)])
    expect_equal(force_intcovar(X[,3], X[,2]), X[-(2:5),c(3,2)])

    # try some permutations
    expect_equal(force_intcovar(X, X[sample(1:nrow(X)),]), X[-(1:5),])
    expect_equal(force_intcovar(X[,-1], X[sample(1:nrow(X)),1]), X[-(1:5),c(2,3,1)])

})

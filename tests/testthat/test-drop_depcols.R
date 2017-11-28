context("drop linearly dependent columns")

test_that("find_lin_indep_cols works", {

    set.seed(20151130)
    x <- cbind(1, sample(0:1, 200, repl=TRUE))

    expect_equal(sort(find_lin_indep_cols(x)), c(1,2))
    expect_equal(sort(find_lin_indep_cols(cbind(x, 1))), c(1,2))
    expect_equal(sort(find_lin_indep_cols(cbind(x, x[,1] + 0.5*x[,2]))), c(2,3))

    X <- matrix(rnorm(1000*10), ncol=10)
    expect_equal(sort(find_lin_indep_cols(X)), 1:10)
    expect_equal(length(find_lin_indep_cols(cbind(rowSums(X), X))), 10)

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
    expect_equal(ncol(drop_depcols(X)), 3)

})

test_that("drop_xcovar works", {

    set.seed(20151202)
    n <- 100
    Xcovar <- sample(0:1, n, replace=TRUE)
    int <- rep(1, n)
    names(Xcovar) <- names(int) <- paste(1:n)

    expect_equal( drop_xcovar(NULL, NULL), NULL)
    expect_equal( drop_xcovar(NULL, Xcovar), Xcovar)
    expect_equal( drop_xcovar(int, NULL), NULL)
    expect_equal( drop_xcovar(int, Xcovar), as.matrix(Xcovar))
    expect_equal( drop_xcovar(Xcovar, Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar), NULL)

    pgm <- sample(0:1, n, replace=TRUE)
    names(pgm) <- paste(1:n)
    Xcovar <- cbind(Xcovar, pgm)

    expect_equal( drop_xcovar(NULL, Xcovar), Xcovar)
    expect_equal( drop_xcovar(int, Xcovar), Xcovar)
    expect_equal( drop_xcovar(Xcovar, Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar[,1]), Xcovar), Xcovar[,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(int, Xcovar[,2]), Xcovar), Xcovar[,1,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar[,1], int), Xcovar), Xcovar[,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar[,2], int), Xcovar), Xcovar[,1,drop=FALSE])

})

test_that("drop_xcovar works with NAs", {

    set.seed(20151202)
    n <- 100
    Xcovar <- sample(0:1, n, replace=TRUE)
    int <- rep(1, n)
    names(Xcovar) <- names(int) <- paste(1:n)
    Xcovar[3:5] <- NA

    expect_equal( drop_xcovar(NULL, Xcovar), Xcovar)
    expect_equal( drop_xcovar(int, Xcovar), as.matrix(Xcovar))
    expect_equal( drop_xcovar(Xcovar, Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar), NULL)

    # try some permutations
    o <- sample(1:n)
    expect_equal( drop_xcovar(int, Xcovar[o]), as.matrix(Xcovar[o]))
    expect_equal( drop_xcovar(Xcovar, Xcovar[o]), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar[o]), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar[o]), NULL)

    # multiple X covar columns
    pgm <- sample(0:1, n, replace=TRUE)
    names(pgm) <- paste(1:n)
    Xcovar <- cbind(Xcovar, pgm)

    expect_equal( drop_xcovar(NULL, Xcovar), Xcovar)
    expect_equal( drop_xcovar(int, Xcovar), Xcovar)
    expect_equal( drop_xcovar(Xcovar, Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar[,1]), Xcovar), Xcovar[,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(int, Xcovar[,2]), Xcovar), Xcovar[,1,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar[,1], int), Xcovar), Xcovar[,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar[,2], int), Xcovar), Xcovar[,1,drop=FALSE])

    # try some permutations
    o <- sample(1:n)
    expect_equal( drop_xcovar(int, Xcovar[o,]), Xcovar[o,])
    expect_equal( drop_xcovar(Xcovar, Xcovar[o,]), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar), Xcovar[o,]), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar[o,]), NULL)
    expect_equal( drop_xcovar(cbind(int, Xcovar[,1]), Xcovar[o,]), Xcovar[o,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(int, Xcovar[,2]), Xcovar[o,]), Xcovar[o,1,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar, int), Xcovar[o,]), NULL)
    expect_equal( drop_xcovar(cbind(Xcovar[,1], int), Xcovar[o,]), Xcovar[o,2,drop=FALSE])
    expect_equal( drop_xcovar(cbind(Xcovar[,2], int), Xcovar[o,]), Xcovar[o,1,drop=FALSE])

})

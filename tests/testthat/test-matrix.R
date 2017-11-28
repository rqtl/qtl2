context("c++ matrix utilities")

test_that("formX_intcovar works", {

    set.seed(20151130)

    # create a prob matrix
    library(qtl)
    data(fake.f2)
    prob <- aperm(calc.genoprob(fake.f2["2",], step=5)$geno[["2"]]$prob[,,-1],
                  c(1,3,2))
    dimnames(prob) <- NULL

    n <- nind(fake.f2)
    addcovar <- cbind(rep(1,n), sample(0:1, n, replace=TRUE))
    intcovar <- addcovar[,2,drop=FALSE]

    X <- formX_intcovar(prob, addcovar, intcovar, 0)
    expected <- cbind(addcovar, prob[,,1], prob[,,1]*intcovar[,1]) # the [,1] makes intcovar an ordinary vector
    expect_equal(X, expected)

    X <- formX_intcovar(prob, addcovar, intcovar, 5)
    expected <- cbind(addcovar, prob[,,6], prob[,,6]*intcovar[,1]) # the [,1] makes intcovar an ordinary vector
    expect_equal(X, expected)

    # no interactive covariates
    expect_equal(formX_intcovar(prob, addcovar, matrix(ncol=0, nrow=n), 2), cbind(addcovar, prob[,,3]))

    # neither interactive nor additive covariates
    expect_equal(formX_intcovar(prob, matrix(ncol=0, nrow=n), matrix(ncol=0, nrow=n), 3), prob[,,4])

    # mismatch in rows
    expect_error(formX_intcovar(prob, addcovar[-n,], intcovar, 0))
    expect_error(formX_intcovar(prob, addcovar, intcovar[-1,], 2))

    # two interactive covariates
    addcovar <- cbind(addcovar, sample(0:1, n, replace=TRUE))

    intcovar <- addcovar[,-1]

    X <- formX_intcovar(prob, addcovar, intcovar, 4)
    expected <- cbind(addcovar, prob[,,5], prob[,,5]*intcovar[,1], prob[,,5]*intcovar[,2])
    expect_equal(X, expected)

})

test_that("expand_genoprobs_intcovar works", {

    set.seed(20151201)
    library(qtl)
    data(fake.f2)
    pr <- aperm(calc.genoprob(fake.f2["2",], step=5)$geno[["2"]]$prob[,,-1],
                c(1,3,2))
    dimnames(pr) <- NULL
    rownames(pr) <- paste(1:nind(fake.f2))

    intcovar <- cbind(sample(0:1, nrow(pr), replace=TRUE),
                      sample(0:1, nrow(pr), replace=TRUE))

    result <- expand_genoprobs_intcovar(pr, intcovar)

    expect_equal(dim(result), dim(pr)*c(1,ncol(intcovar)+1,1))
    expect_equal(result[,1:2,] , pr)
    expect_equal(result[,3:4,], pr*intcovar[,1])
    expect_equal(result[,5:6,], pr*intcovar[,2])

    # single interactive covariate
    result2 <- expand_genoprobs_intcovar(pr, intcovar[,1,drop=FALSE])

    expect_equal(dim(result2), dim(pr)*c(1,2,1))
    expect_equal(result2[,1:2,] , pr)
    expect_equal(result2[,3:4,], pr*intcovar[,1])

})


test_that("weighted_matrix works", {

    set.seed(20151201)
    n <- 100
    p <- 10
    X <- matrix(rnorm(n*p), ncol=p)
    w <- runif(n, 1, 4)

    result <- weighted_matrix(X, w)

    expect_equal(result, X*w)
    for(i in 1:p) expect_equal(result[,i], X[,i]*w)

})

test_that("weighted_3darray works", {

    set.seed(20151201)
    n <- 100
    p <- 3
    q <- 8
    X <- array(rnorm(n*p*q), dim=c(n, p, q))
    w <- runif(n, 1, 4)

    result <- weighted_3darray(X, w)

    expect_equal(result, X*w)
    for(i in 1:p)
        for(j in 1:q)
            expect_equal(result[,i,j], X[,i,j]*w)

})

test_that("matrix_x_matrix works", {
    X <- matrix(rnorm(200), ncol=10)
    Y <- matrix(rnorm(50), nrow=10)
    expect_equal(matrix_x_matrix(X, Y), X %*% Y)
})

test_that("matrix_x_vector works", {
    X <- matrix(rnorm(200), ncol=10)
    Y <- matrix(rnorm(30), nrow=10)
    expect_equal(matrix_x_vector(X, Y[,1]), (X %*% Y[,1])[,1])
    expect_equal(matrix_x_vector(X, Y[,2]), (X %*% Y[,2])[,1])
    expect_equal(matrix_x_vector(X, Y[,3]), (X %*% Y[,3])[,1])
})

test_that("matrix_x_3darray works", {
    p <- 25
    n <- 50
    g <- 8
    pos <- 100
    X <- matrix(rnorm(p*n), ncol=n)
    Y <- array(rnorm(n*g*pos), dim=c(n, g, pos))

    expect <- array(dim=c(nrow(X), dim(Y)[2:3]))
    for(i in 1:dim(expect)[3])
        expect[,,i] <- X %*% Y[,,i]

    expect_equal(matrix_x_3darray(X, Y), expect)

})

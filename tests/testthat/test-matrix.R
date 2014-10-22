context("c++ matrix utilities")

test_that("cbind_imatrix works", {

    set.seed(54028069)
    x1 <- matrix(sample(1:10000, 100), nrow=10)
    x2 <- matrix(sample(1:10000, 200), nrow=10)
    x3 <- matrix(sample(1:10000, 50), nrow=10)

    expect_equal(cbind_imatrix(x1, x2), cbind(x1, x2))
    expect_equal(cbind_3imatrix(x1, x2, x3), cbind(x1, x2, x3))

    expect_error(cbind_imatrix(x1, t(x2)))
    expect_error(cbind_3imatrix(x1, t(x2), x3))
    expect_error(cbind_3imatrix(x1, x2, t(x3)))

})

test_that("cbind_nmatrix works", {

    set.seed(54028069)
    x1 <- matrix(runif(100, 0, 100), nrow=10)
    x2 <- matrix(runif(200, 0, 100), nrow=10)
    x3 <- matrix(runif( 50, 0, 100), nrow=10)

    expect_equal(cbind_nmatrix(x1, x2), cbind(x1, x2))
    expect_equal(cbind_3nmatrix(x1, x2, x3), cbind(x1, x2, x3))

    expect_error(cbind_nmatrix(x1, t(x2)))
    expect_error(cbind_3nmatrix(x1, t(x2), x3))
    expect_error(cbind_3nmatrix(x1, x2, t(x3)))

})

test_that("rbind_imatrix works", {

    set.seed(54028069)
    x1 <- matrix(sample(1:10000, 100), ncol=10)
    x2 <- matrix(sample(1:10000, 200), ncol=10)
    x3 <- matrix(sample(1:10000, 50), ncol=10)

    expect_equal(rbind_imatrix(x1, x2), rbind(x1, x2))
    expect_equal(rbind_3imatrix(x1, x2, x3), rbind(x1, x2, x3))

    expect_error(rbind_imatrix(x1, t(x2)))
    expect_error(rbind_3imatrix(x1, t(x2), x3))
    expect_error(rbind_3imatrix(x1, x2, t(x3)))

})

test_that("rbind_nmatrix works", {

    set.seed(54028069)
    x1 <- matrix(runif(100, 0, 100), ncol=10)
    x2 <- matrix(runif(200, 0, 100), ncol=10)
    x3 <- matrix(runif( 50, 0, 100), ncol=10)

    expect_equal(rbind_nmatrix(x1, x2), rbind(x1, x2))
    expect_equal(rbind_3nmatrix(x1, x2, x3), rbind(x1, x2, x3))

    expect_error(rbind_nmatrix(x1, t(x2)))
    expect_error(rbind_3nmatrix(x1, t(x2), x3))
    expect_error(rbind_3nmatrix(x1, x2, t(x3)))

})

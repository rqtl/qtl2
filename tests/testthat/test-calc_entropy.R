context("Calculate entropy")

test_that("calc_entropy works in a simple case", {

    p <- list("1"=array(c(1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, 1,
                          0.25, 0.25, 0.25, 0.25,
                          0.5, 0.5, 0, 0,
                          0.5, 0, 0.5, 0,
                          0.5, 0, 0, 0.5,
                          0, 0.5, 0.5, 0,
                          0, 0.5, 0, 0.5,
                          0, 0, 0.5, 0.5,
                          0.25, 0.25, 0.5, 0), dim=c(4, 2, 6)))
    p[[1]] <- aperm(p[[1]], c(2,1,3))
    dimnames(p[[1]]) <- list(c("ind1", "ind2"), letters[1:4], LETTERS[1:6])

    expected <- list("1"=cbind(c(0,0),
                                c(0,0),
                                c(2,1),
                                c(1,1),
                                c(1,1),
                                c(1,1.5)))
    dimnames(expected[[1]]) <- dimnames(p[[1]])[c(1,3)]

    expect_equal(calc_entropy(p), expected)

})

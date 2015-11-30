context("find matching cols")

test_that("find_matchin_cols works", {

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

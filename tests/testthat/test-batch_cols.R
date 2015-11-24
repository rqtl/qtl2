context("batch columns by NAs")

test_that("batch_cols works", {

    x <- rbind(c( 1,  2,  3),
               c( 4,  5,  6),
               c( 7, NA,  8),
               c(NA, NA,  9),
               c(10, 11, 12))

    expect_equal(batch_cols(x),
                 list(list(cols=1, keep=c(TRUE, TRUE, TRUE, FALSE, TRUE)),
                      list(cols=2, keep=c(TRUE, TRUE, FALSE, FALSE, TRUE)),
                      list(cols=3, keep=c(TRUE, TRUE, TRUE, TRUE, TRUE)))
                 )

    expect_equal(batch_cols( rbind(c(1,2,3)) ),
                 list(list(cols=c(1,2,3), keep=TRUE))
                 )

    expect_equal(batch_cols( rbind(1, 2, NA, 3) ),
                 list(list(cols=1, keep=c(TRUE, TRUE, FALSE, TRUE)))
                 )

    y <- rbind(c( 1,  2,  3, 13, 16),
               c( 4,  5,  6, 14, 17),
               c( 7, NA,  8, NA, 18),
               c(NA, NA, NA, NA, 19),
               c(10, 11, 12, 15, 20))
    expect_equal(batch_cols(y),
                 list(list(cols=c(1,3), keep=c(TRUE, TRUE, TRUE, FALSE, TRUE)),
                      list(cols=c(2,4), keep=c(TRUE, TRUE, FALSE, FALSE, TRUE)),
                      list(cols=5, keep=rep(TRUE, 5)))
                 )
})

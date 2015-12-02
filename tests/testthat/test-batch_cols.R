context("batch columns by NAs")

test_that("batch_cols works", {

    x <- rbind(c( 1,  2,  3),
               c( 4,  5,  6),
               c( 7, NA,  8),
               c(NA, NA,  9),
               c(10, 11, 12))

    expect_equal(batch_cols(x),
                 list(list(cols=3, omit=numeric(0)),
                      list(cols=1, omit=4),
                      list(cols=2, omit=c(3,4)))
                 )

    expect_equal(batch_cols( rbind(c(1,2,3)) ),
                 list(list(cols=c(1,2,3), omit=numeric(0)))
                 )

    expect_equal(batch_cols( rbind(1, 2, NA, 3) ),
                 list(list(cols=1, omit=3))
                 )

    y <- rbind(c( 1,  2,  3, 13, 16),
               c( 4,  5,  6, 14, 17),
               c( 7, NA,  8, NA, 18),
               c(NA, NA, NA, NA, 19),
               c(10, 11, 12, 15, 20))
    expect_equal(batch_cols(y),
                 list(list(cols=5, omit=numeric(0)),
                      list(cols=c(1,3), omit=4),
                      list(cols=c(2,4), omit=c(3,4)))

                 )
})

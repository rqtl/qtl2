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

test_that("batch_vec works", {

    # one big batch
    expect_equal(batch_vec(1:100, 100), list(1:100))

    # one element
    expect_equal(batch_vec(1, 1), list(1))
    expect_equal(batch_vec(1, 8), list(1))

    # fewer elements than batch size
    expect_equal(batch_vec(1:5, 8), list(1:5))
    expect_equal(batch_vec(1:30, 32), list(1:30))

    # equal elements in each core
    expect <- as.list(as.data.frame(matrix(1:10, nrow=5, byrow=FALSE)))
    names(expect) <- NULL
    expect_equal(batch_vec(1:10, 5), expect)

    expect <- as.list(as.data.frame(matrix(1:10, nrow=2, byrow=FALSE)))
    names(expect) <- NULL
    expect_equal(batch_vec(1:10, 2), expect)

    expect <- as.list(as.data.frame(matrix(1:144, nrow=18, byrow=FALSE)))
    names(expect) <- NULL
    expect_equal(batch_vec(1:144, 18), expect)

    expect <- as.list(as.data.frame(matrix(1:144, nrow=8, byrow=FALSE)))
    names(expect) <- NULL
    expect_equal(batch_vec(1:144, 8), expect)

    # non-equal
    expect_equal(batch_vec(101:110, 3),
                      list(101:103, 104:106, 107:108, 109:110))
    expect_equal(batch_vec(1:35, 6),
                      list(1:6, 7:12, 13:18, 19:24, 25:30, 31:35))
    expect_equal(batch_vec((1:16)+3, 9),
                      list((1:8)+3, (9:16)+3))

    set.seed(20151202)
    vec <- sample((1:17)+4)
    expect_equal(batch_vec(vec, 9),
                      list(vec[1:9], vec[10:17]))

})


test_that("batch_cols with max_batch works", {

    set.seed(20151202)

    n <- 100
    p <- 1000
    x <- matrix(rnorm(n*p), ncol=p)
    expect_equal(batch_cols(x),
                 list(list(cols=1:p, omit=numeric(0))))

    expect <- vector("list", 4)
    for(i in 1:4)
        expect[[i]] <- list(cols=1:250+(i-1)*250, omit=numeric(0))
    expect_equal(batch_cols(x, max_batch=250), expect)

    x[16:20, 8:12] <- NA
    expect <- vector("list", 5)
    expect[[1]] <- list(cols=8:12, omit=16:20)
    rest <- (1:1000)[-(8:12)]
    expect[[2]] <- list(cols=rest[1:249], omit=numeric(0))
    expect[[3]] <- list(cols=rest[249+1:249], omit=numeric(0))
    expect[[4]] <- list(cols=rest[2*249+1:249], omit=numeric(0))
    expect[[5]] <- list(cols=rest[3*249+1:248], omit=numeric(0))
    expect_equal(batch_cols(x, max_batch=250), expect)

})

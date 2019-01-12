context("cbind_expand")

test_that("cbind_expand works", {

    set.seed(20151109)
    x <- data.frame(x=sample(LETTERS[1:5]), y=sample(letters[6:10]), stringsAsFactors=FALSE)
    rownames(x) <- c("1", "3", "4", "5", "6")
    y <- data.frame(a=runif(4), b=rnorm(4))
    rownames(y) <- c("5", "2", "9", "4")

    expected <- data.frame(x=c(x$x, NA, NA), y=c(x$y, NA, NA),
                           a=rep(NA, 7), b=rep(NA, 7),
                           stringsAsFactors=FALSE)
    rownames(expected) <- c("1", "3", "4", "5", "6", "2", "9")
    expected[rownames(y),"a"] <- y$a
    expected[rownames(y),"b"] <- y$b


    expect_equal(cbind_expand(x, y), expected)

    # reverse inputs
    expect_equal(cbind_expand(y, x), expected[c(4,6,7,3,1,2,5), c(3,4,1,2)])

    # duplicate inputs
    expect_equal(cbind_expand(x, x), cbind(x, x))
    expect_equal(cbind_expand(y, y), cbind(y, y))

    # expect error due to duplicate row names
    expect_error(cbind_expand( rbind("1"=1:3, "1"=4:6), rbind("1"=4:9) ) )


    # add a third input
    z <- data.frame(c=sample(LETTERS[11:17]), d=sample(1L:5L, 7, replace=TRUE),
                    stringsAsFactors=FALSE)
    rownames(z) <- c("3", "1", "9", "4", "6", "5", "8")

    expected <- data.frame(x=c(x$x, NA, NA, NA), y=c(x$y, NA, NA, NA),
                           a=rep(NA, 8), b=rep(NA, 8),
                           c=rep("", 8), d=rep(NA, 8),
                           stringsAsFactors=FALSE)
    rownames(expected) <- c("1", "3", "4", "5", "6", "2", "9", "8")
    expected[rownames(y),"a"] <- y$a
    expected[rownames(y),"b"] <- y$b
    expected[rownames(z),"d"] <- z$d
    expected[expected$c=="", "c"] <- NA
    expected[rownames(z),"c"] <- z$c
    expect_equal(cbind_expand(x, y, z), expected)


    # a few more simple tests
    df1 <- data.frame(x=c(1,2,3,NA,4), y=c(5,8,9,10,11), row.names=c("A", "B", "C", "D", "E"))
    df2 <- data.frame(w=c(7,8,0,9,10), z=c(6,NA,NA,9,10), row.names=c("A", "B", "F", "C", "D"))

    df1n2 <- data.frame(x=c(1,2,3,NA,4,NA), y=c(5,8,9,10,11,NA),w=c(7,8,9,10,NA,0),z=c(6,NA,9,10,NA,NA),
                        row.names=c("A", "B", "C", "D", "E","F"))
    df2n1 <- df1n2[c(1,2,6,3:5), c(3:4,1:2)]

    expect_equal(cbind_expand(df1, df2), df1n2)
    expect_equal(cbind_expand(df2, df1), df2n1)

    expect_equal(cbind_expand(df1, df1), cbind(df1, df1))
    expect_equal(cbind_expand(df2, df2), cbind(df2, df2))

})

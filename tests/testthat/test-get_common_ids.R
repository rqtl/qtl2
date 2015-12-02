context("get common ids from set of objects")

test_that("get_common_ids works", {

    x <- matrix(1:20, nrow=5, ncol=4)
    rownames(x) <- LETTERS[1:5]

    y <- 1:5
    names(y) <- c("A", "C", "D", "E", "F")

    z <- array(1:(5*2*3), dim=c(5, 2, 3))
    dimnames(z) <- list(c("G", "C", "D", "G", "A"),
                        c("col1", "col2"),
                        c("depth1", "depth2", "depth3"))

    # error due to duplicate name
    expect_error(get_common_ids(x, y, z))
    expect_error(get_common_ids(z, y, x))
    expect_error(get_common_ids(x, z, y))

    rownames(z)[4] <- "H" # fix dup
    expect_equal(get_common_ids(x, y, z), c("A", "C", "D"))
    expect_equal(get_common_ids(x, z, y), c("A", "C", "D"))
    expect_equal(get_common_ids(z, y, x), c("C", "D", "A"))

    rownames(z)[4] <- "H" # fix dup
    chvec <- c("A", "I", "D", "E")
    expect_equal(get_common_ids(x, y, z, chvec), c("A", "D"))
    expect_equal(get_common_ids(x, z, chvec, y), c("A", "D"))
    expect_equal(get_common_ids(z, chvec, y, x), c("D", "A"))

    # try data.frame
    x <- as.data.frame(x)
    rownames(z)[4] <- "G"
    expect_error(get_common_ids(x, y, z)) # duplicate names
    expect_error(get_common_ids(z, y, x)) # duplicate names
    expect_error(get_common_ids(x, z, y)) # duplicate names

    rownames(z)[4] <- "H"
    expect_equal(get_common_ids(x, y, z), c("A", "C", "D"))
    expect_equal(get_common_ids(x, z, y), c("A", "C", "D"))
    expect_equal(get_common_ids(z, y, x), c("C", "D", "A"))

})

test_that("get_common_ids works when complete.cases=TRUE", {

    n <- 10
    x <- matrix(rnorm(10*n), ncol=10)
    rownames(x) <- LETTERS[1:10]

    y <- matrix(rnorm(5*(n-1)), ncol=5)
    rownames(y) <- LETTERS[(1:(n-1))+2]

    z <- LETTERS[c(3:8, 15:18)]

    expected <- LETTERS[3:8]
    expect_equal(get_common_ids(x, y, z, complete.cases=TRUE), expected)
    expect_equal(get_common_ids(y, x, z, complete.cases=TRUE), expected)
    expect_equal(get_common_ids(z, y, x, complete.cases=TRUE), expected)

    # add some missing values
    x[1,c(5,7)] <- NA    # A
    x[3,c(2,8,9)] <- NA  # C
    x[5,4] <- NA         # E
    y[3,5] <- NA         # C
    y[6,] <- NA          # H

    expected <- c("D", "F", "G")
    expect_equal(get_common_ids(x, y, z, complete.cases=TRUE), expected)
    expect_equal(get_common_ids(y, x, z, complete.cases=TRUE), expected)
    expect_equal(get_common_ids(z, y, x, complete.cases=TRUE), expected)

})

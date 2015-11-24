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

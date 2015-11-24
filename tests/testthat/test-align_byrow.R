context("align objects by row names")

test_that("align_byrow works", {

    x <- matrix(1:20, nrow=5, ncol=4)
    rownames(x) <- LETTERS[1:5]

    y <- 1:5
    names(y) <- c("A", "C", "D", "E", "F")

    z <- array(1:(5*2*3), dim=c(5, 2, 3))
    dimnames(z) <- list(c("G", "C", "D", "G", "A"),
                        c("col1", "col2"),
                        c("depth1", "depth2", "depth3"))

    # error due to duplicate name
    expect_error(align_byrow(x, y, z))
    expect_error(align_byrow(z, y, x))
    expect_error(align_byrow(x, z, y))
    expect_error(align_byrow(x, y, z, return_index=TRUE))
    expect_error(align_byrow(z, y, x, return_index=TRUE))
    expect_error(align_byrow(x, z, y, return_index=TRUE))

    rownames(z)[4] <- "H" # fix dup
    expect_equal(align_byrow(x, y, z, return_index=TRUE), c("A", "C", "D"))
    expect_equal(align_byrow(x, z, y, return_index=TRUE), c("A", "C", "D"))
    expect_equal(align_byrow(z, y, x, return_index=TRUE), c("C", "D", "A"))

    id <- c("A", "C", "D")
    expect_equal(align_byrow(x, y, z), list(x[id,,drop=FALSE], y[id], z[id,,,drop=FALSE]))
    expect_equal(align_byrow(x, z, y), list(x[id,,drop=FALSE], z[id,,,drop=FALSE], y[id]))
    id <- c("C", "D", "A")
    expect_equal(align_byrow(z, y, x), list(z[id,,,drop=FALSE], y[id], x[id,,drop=FALSE]))

    # try data.frame
    x <- as.data.frame(x)
    rownames(z)[4] <- "G"
    expect_error(align_byrow(x, y, z)) # duplicate names
    expect_error(align_byrow(z, y, x)) # duplicate names
    expect_error(align_byrow(x, z, y)) # duplicate names
    expect_error(align_byrow(x, y, z, return_index=TRUE)) # duplicate names
    expect_error(align_byrow(z, y, x, return_index=TRUE)) # duplicate names
    expect_error(align_byrow(x, z, y, return_index=TRUE)) # duplicate names

    rownames(z)[4] <- "H"
    expect_equal(align_byrow(x, y, z, return_index=TRUE), c("A", "C", "D"))
    expect_equal(align_byrow(x, z, y, return_index=TRUE), c("A", "C", "D"))
    expect_equal(align_byrow(z, y, x, return_index=TRUE), c("C", "D", "A"))

    id <- c("A", "C", "D")
    expect_equal(align_byrow(x, y, z), list(x[id,,drop=FALSE], y[id], z[id,,,drop=FALSE]))
    expect_equal(align_byrow(x, z, y), list(x[id,,drop=FALSE], z[id,,,drop=FALSE], y[id]))
    id <- c("C", "D", "A")
    expect_equal(align_byrow(z, y, x), list(z[id,,,drop=FALSE], y[id], x[id,,drop=FALSE]))

})

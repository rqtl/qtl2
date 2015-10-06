context("find_common_ids")

test_that("find_common_ids works", {

    set.seed(20151006)
    first <- LETTERS[1:10]
    second <- sample(LETTERS[c(1:4,7:12)])

    m12 <- match(first, second)
    expected <- list(first=c(1:4, 7:10),
                     second=m12[!is.na(m12)],
                     firstonly=5:6,
                     secondonly=which(second %in% LETTERS[11:12]))
    names(expected$first) <- first[expected$first]
    names(expected$second) <- second[expected$second]
    names(expected$firstonly) <- first[expected$firstonly]
    names(expected$secondonly) <- second[expected$secondonly]

    expect_equal(find_common_ids(first, second), expected)
    expect_equal(first[expected$first], second[expected$second])

})

test_that("find_common_ids works with no overlap", {

    set.seed(20151006)
    first <- LETTERS[1:10]
    second <- LETTERS[11:20]

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=1:10)
    names(expected$firstonly) <- first
    names(expected$secondonly) <- second

    expect_equal(find_common_ids(first, second), expected)

})

test_that("find_common_ids works with null vectors", {

    set.seed(20151006)
    first <- LETTERS[1:10]
    second <- character(0)

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=numeric(0))
    names(expected$firstonly) <- first

    expect_equal(find_common_ids(first, second), expected)

    # switch the role of the two
    expected$secondonly <- expected$firstonly
    expected$firstonly <- numeric(0)
    expect_equal(find_common_ids(second, first), expected)

    # same, when second is NULL
    second <- NULL

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=numeric(0))
    names(expected$firstonly) <- first

    expect_equal(find_common_ids(first, second), expected)

    # switch the role of the two (when second is NULL)
    expected$secondonly <- expected$firstonly
    expected$firstonly <- numeric(0)
    expect_equal(find_common_ids(second, first), expected)

    # both NULL
    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=numeric(0),
                     secondonly=numeric(0))
    expect_equal(find_common_ids(NULL, NULL), expected)

})

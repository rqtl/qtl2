context("find_common_ids and match_ids")

test_that("match_ids and find_common_ids work", {

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

    expect_equal(match_ids(first, second), expected)
    expect_equal(first[expected$first], second[expected$second])

    expect_equal(find_common_ids(first, second), first[expected$first])

})

test_that("match_ids and find_common_ids work with no overlap", {

    set.seed(20151006)
    first <- LETTERS[1:10]
    second <- LETTERS[11:20]

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=1:10)
    names(expected$firstonly) <- first
    names(expected$secondonly) <- second

    expect_equal(match_ids(first, second), expected)
    expect_equal(first[expected$first], second[expected$second])

    expect_equal(find_common_ids(first, second), first[expected$first])

})

test_that("match_ids and find_common_ids work with null vectors", {

    set.seed(20151006)
    first <- LETTERS[1:10]
    second <- character(0)

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=numeric(0))
    names(expected$firstonly) <- first

    expect_equal(match_ids(first, second), expected)
    expect_equal(first[expected$first], second[expected$second])
    expect_equal(find_common_ids(first, second), first[expected$first])

    # switch the role of the two
    expected$secondonly <- expected$firstonly
    expected$firstonly <- numeric(0)
    expect_equal(match_ids(second, first), expected)
    expect_equal(second[expected$first], first[expected$second])
    expect_equal(find_common_ids(second, first), second[expected$first])

    # same, when second is NULL
    second <- NULL

    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=1:10,
                     secondonly=numeric(0))
    names(expected$firstonly) <- first

    expect_equal(match_ids(first, second), expected)
    expect_equal(find_common_ids(first, second), first[expected$first])

    # switch the role of the two (when second is NULL)
    expected$secondonly <- expected$firstonly
    expected$firstonly <- numeric(0)
    expect_equal(match_ids(second, first), expected)

    # both NULL
    expected <- list(first=numeric(0),
                     second=numeric(0),
                     firstonly=numeric(0),
                     secondonly=numeric(0))
    expect_equal(match_ids(NULL, NULL), expected)
    expect_equal(find_common_ids(NULL, NULL), character(0))

    expect_equal(match_ids(character(0), character(0)), expected)
    expect_equal(find_common_ids(character(0), character(0)), character(0))

})


test_that("find_common_ids with <2 inputs", {

    expect_error(find_common_ids())

    expect_equal(find_common_ids(LETTERS[1:5]), LETTERS[1:5])
    expect_equal(find_common_ids(NULL), character(0))
    expect_equal(find_common_ids(character(0)), character(0))

})


test_that("find_common_ids with >2 inputs", {

    set.seed(20151008)
    first <- LETTERS[1:10]
    second <- sample(LETTERS[3:12])
    third <- sample(LETTERS[5:14])

    expect_equal(find_common_ids(first, second, third), LETTERS[5:10])
    expect_equal(sort(find_common_ids(second, first, third)), LETTERS[5:10])
    expect_equal(sort(find_common_ids(third, second, first)), LETTERS[5:10])

    expect_equal(find_common_ids(third, second, NULL, first), character(0))

    fourth <- sample(c("H", "E", "L", "P"))
    expect_equal(sort(find_common_ids(third, second, fourth, first)), c("E", "H"))

})

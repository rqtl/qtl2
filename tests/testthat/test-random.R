
context("random/permutation")

test_that("random_int looks correct", {

    set.seed(71516331)
    n <- 10000
    lo <- 5
    hi <- 21
    x <- floor(runif(n, lo, hi+1))
    storage.mode(x) <- "integer"
    set.seed(71516331)
    y <- random_int(n, lo, hi)
    expect_equal(x, y)
})

test_that("permutations look correct", {

    set.seed(701777)
    n.runs <- 10
    size <- 100
    x <- runif(size)
    # random permutations of x
    result1 <- replicate(n.runs, permute_nvector(x))
    # same stuff after sorting?
    sorted1 <- apply(result1, 2, sort)
    sortedx <- sort(x)
    for(i in 1:n.runs)
        expect_equal(sorted1[,i], sortedx)

    set.seed(40451392)
    x <- 1:size
    # random permutations of x
    result2 <- replicate(n.runs, permute_ivector(x))
    # same stuff after sorting?
    sorted2 <- apply(result2, 2, sort)
    for(i in 1:n.runs)
        expect_equal(sorted2[,i], x)

    set.seed(40451392)
    # random permutations of 0:(size-1)
    result3 <- replicate(n.runs, get_permutation(size))
    # same stuff as result above?
    for(i in 1:n.runs)
        expect_equal(result2[,i], result3[,i]+1)

    set.seed(701777)
    x <- runif(size)
    result4 <- get_permutations(n.runs, x)
    for(i in 1:n.runs)
        expect_equal(result1[,i], result4[,i])

})

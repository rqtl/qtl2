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
    n_runs <- 10
    size <- 100
    x <- runif(size)
    # random permutations of x
    result1 <- permute_nvector(n_runs, x)
    # same stuff after sorting?
    sorted1 <- apply(result1, 2, sort)
    sortedx <- sort(x)
    for(i in 1:n_runs)
        expect_equal(sorted1[,i], sortedx)

    set.seed(40451392)
    x <- 1:size
    # random permutations of x
    result2 <- permute_ivector(n_runs, x)
    # same stuff after sorting?
    sorted2 <- apply(result2, 2, sort)
    for(i in 1:n_runs)
        expect_equal(sorted2[,i], x)

    set.seed(40451392)
    # random permutations of 0:(size-1)
    result3 <- replicate(n_runs, get_permutation(size))
    # same stuff as result above?
    for(i in 1:n_runs)
        expect_equal(result2[,i], result3[,i]+1)

})

test_that("stratified permutations for integers works", {

    vals <- sample(1:10000, 250)
    strata <- sample(0:4, 250, replace=TRUE)

    set.seed(22632617)
    z <- permute_ivector_stratified(100, vals, strata)

    expect_true( all(apply(z, 2, function(a,b) all(sort(a)==b), sort(vals))) )

    for(i in unique(strata))
        expect_true( all(apply(z[strata==i,], 2, function(a,b) all(sort(a)==b), sort(vals[strata==i]))) )

    set.seed(22632617)
    z2 <- permute_ivector_stratified(100, vals, strata, max(strata)+1)

    expect_equal(z, z2)

    # error: n_strata is too small
    expect_error(permute_ivector_stratified(100, vals, strata, 3),
                 'strata should be in \\[0, n_strata)')

    # error: vals and strata are different lengths
    expect_error(permute_ivector_stratified(100, vals[1:200], strata))

})

test_that("stratified permutations for dobules works", {

    vals <- runif(250, 0, 100)
    strata <- sample(0:4, 250, replace=TRUE)

    set.seed(22632617)
    z <- permute_nvector_stratified(100, vals, strata)

    expect_true( all(apply(z, 2, function(a,b) all(sort(a)==b), sort(vals))) )

    for(i in unique(strata))
        expect_true( all(apply(z[strata==i,], 2, function(a,b) all(sort(a)==b), sort(vals[strata==i]))) )

    set.seed(22632617)
    z2 <- permute_nvector_stratified(100, vals, strata, max(strata)+1)

    expect_equal(z, z2)

    # error: n_strata is too small
    expect_error(permute_nvector_stratified(100, vals, strata, 3),
                 'strata should be in \\[0, n_strata)')

    # error: vals and strata are different lengths
    expect_error(permute_nvector_stratified(100, vals[1:200], strata))
})

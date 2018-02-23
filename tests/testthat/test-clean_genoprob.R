context("Clean genotype probabilities")

test_that("clean_genoprob works", {

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c("19", "X")]

    map <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, map, error_prob=0.002, map_function="c-f")

    # there are no *really* small values here
    pr_clean <- clean_genoprob(pr)
    expect_equal(pr_clean, pr)

    # use a bigger tolerance for the individual values
    pr_clean <- clean_genoprob(pr, 0.001, 0.01)

    for(i in seq_along(pr_clean)) {
        # ensure that things still sum to 1
        expect_true( max(abs(apply(pr_clean[[i]], c(1,3), sum) - 1)) < 1e-12 )
    }

    # create a column that is largely missing
    pr[[1]][,3,20] <- runif(nrow(pr[[1]]), 0, 1e-4)
    # re-scale so values sum to 1
    pr[[1]][,,20] <- pr[[1]][,,20] / rowSums(pr[[1]][,,20])

    pr_clean <- clean_genoprob(pr, column_threshold=0.01)
    expect_true( all(pr_clean[[1]][,3,20] == 0) )
    expect_equal( pr_clean[[1]][,,-20], pr[[1]][,,-20] )
    expect_equal( pr_clean[[2]], pr[[2]] )

    expect_true( max(abs(apply(pr_clean[[1]], c(1,3), sum) - 1)) < 1e-12 )

    # try this with a subset of individuals
    # create another column that is largely missing, for a subset of individuals
    set.seed(20180223)
    ind <- sample(ind_ids(iron), 50)
    pr[[1]][ind,1,15] <- runif(length(ind), 0, 1e-4)
    # re-scale so values sum to 1
    pr[[1]][,,15] <- pr[[1]][,,15] / rowSums(pr[[1]][,,15])


    pr_clean <- clean_genoprob(pr, column_threshold=0.01, ind=ind)
    other_ind <- ind_ids(iron)[!(ind_ids(iron) %in% ind)]
    expect_equal(pr_clean[other_ind,], pr[other_ind,])
    expect_true( max(abs(apply(pr_clean[[1]], c(1,3), sum) - 1)) < 1e-12 )
    expect_true( all(pr_clean[[1]][ind,1,15] == 0) )

    expect_equal( pr_clean[[1]][,,-c(15,20)], pr[[1]][,,-c(15,20)] )
    expect_equal( pr_clean[[2]], pr[[2]] )
})

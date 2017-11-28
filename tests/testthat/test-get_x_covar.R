context("get X covariates")

test_that("get_x_covar for riself", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    expect_equal(get_x_covar(grav2), NULL)

})

test_that("get_x_covar for intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    n <- n_ind(iron)

    expected <- matrix(0.0, ncol=2, nrow=n)
    expected[!iron$is_female, 1] <- 1.0
    expected[iron$is_female & iron$cross_info==1, 2] <- 1.0
    dimnames(expected) <- list(names(iron$is_female), c("sex", "direction"))

    expect_equal(get_x_covar(iron), expected)

    # all male
    tmp <- iron
    tmp$is_female <- rep(FALSE, length(tmp$is_female))
    expect_equal(get_x_covar(tmp), NULL)

    # all female
    tmp <- iron
    tmp$is_female <- rep(TRUE, length(tmp$is_female))
    expected <- matrix(0.0, ncol=1, nrow=n)
    expected[iron$cross_info==1, 1] <- 1.0
    dimnames(expected) <- list(names(tmp$is_female), "direction")
    expect_equal(get_x_covar(tmp), expected)

    # all female, forward direction
    tmp <- iron
    tmp$is_female <- rep(TRUE, length(tmp$is_female))
    tmp$cross_info[tmp$cross_info==1] <- 0
    expect_equal(get_x_covar(tmp), NULL)

    # all female, reverse direction
    tmp <- iron
    tmp$is_female <- rep(TRUE, length(tmp$is_female))
    tmp$cross_info[tmp$cross_info==0] <- 1
    expect_equal(get_x_covar(tmp), NULL)

    # both sexes, one direction
    tmp <- iron
    tmp$cross_info[tmp$cross_info==1] <- 0
    expected <- matrix(0.0, ncol=1, nrow=n)
    expected[!iron$is_female] <- 1
    dimnames(expected) <- list(names(tmp$is_female), "sex")
    expect_equal(get_x_covar(tmp), expected)

})

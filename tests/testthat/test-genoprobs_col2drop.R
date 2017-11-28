context("genoprobs columns to drop")

test_that("genoprobs_cols2drop works", {

    # load example data from R/qtl
    library(qtl)
    data(fake.f2)
    fake.f2 <- fake.f2[c(18,19,"X"),] # subset chromosomes
    fake.f2 <- convert2cross2(fake.f2)

    # calculate genoprobs
    pr <- calc_genoprob(fake.f2)

    # just X chromosome
    expect_equal(genoprobs_col2drop(pr$X), numeric(0))

    # all chromosomes
    expected <- vector("list", length(pr))
    names(expected) <- names(pr)
    for(i in seq(along=expected)) expected[[i]] <- numeric(0)
    expect_equal(genoprobs_col2drop(pr), expected)
    expect_equal(genoprobs_col2drop(pr, FALSE), expected)

    #####
    # all males
    fake.f2$is_female[fake.f2$is_female] <- FALSE
    pr <- calc_genoprob(fake.f2)

    # just X chromosome
    expect_equal(genoprobs_col2drop(pr$X), 1:4)

    # all chromosomes
    expected$X <- 1:4
    expect_equal(genoprobs_col2drop(pr), expected)
    expect_equal(genoprobs_col2drop(pr, FALSE), expected)

    # test Xonly=FALSE
    expected2 <- expected;expected2$X <- numeric(0)
    prA <- pr
    attr(prA, "is_x_chr") <- rep(FALSE, length(prA))
    expect_equal(genoprobs_col2drop(prA), expected2)
    expect_equal(genoprobs_col2drop(prA, FALSE), expected)

    #####
    # all females
    fake.f2$is_female[!fake.f2$is_female] <- TRUE
    pr <- calc_genoprob(fake.f2)

    # just X chromosome
    expect_equal(genoprobs_col2drop(pr$X), 5:6)

    # all chromosomes
    expected$X <- 5:6
    expect_equal(genoprobs_col2drop(pr), expected)
    expect_equal(genoprobs_col2drop(pr, FALSE), expected)

    # test Xonly=FALSE
    expected2 <- expected;expected2$X <- numeric(0)
    prA <- pr
    attr(prA, "is_x_chr") <- rep(FALSE, length(prA))
    expect_equal(genoprobs_col2drop(prA), expected2)
    expect_equal(genoprobs_col2drop(prA, FALSE), expected)

})

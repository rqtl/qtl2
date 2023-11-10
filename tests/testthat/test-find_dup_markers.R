context("find duplicate markers")

test_that("find_dup_markers matches qtl::findDupMarkers", {

    library(qtl)
    data(hyper)

    hyper2 <- convert2cross2(hyper)

    # exact only, not adjacent only
    set.seed(20231110)
    dup2 <- find_dup_markers(hyper2)

    set.seed(20231110)
    dup <- qtl::findDupMarkers(hyper)

    expect_equal(dup2, dup)

    # exact only, adjacent only
    set.seed(20231110)
    dup2 <- find_dup_markers(hyper2, adjacent_only=TRUE)

    set.seed(20231110)
    dup <- qtl::findDupMarkers(hyper, adjacent.only=TRUE)

    expect_equal(dup2, dup)

    # not exact only, not adjacent only
    set.seed(20231110)
    dup2 <- find_dup_markers(hyper2, exact_only=FALSE)

    set.seed(20231110)
    dup <- qtl::findDupMarkers(hyper, exact.only=FALSE)

    expect_equal(dup2, dup)

    # not exact only, adjacent only
    set.seed(20231110)
    dup2 <- find_dup_markers(hyper2, exact_only=FALSE, adjacent_only=TRUE)

    set.seed(20231110)
    dup <- qtl::findDupMarkers(hyper, exact.only=FALSE, adjacent.only=TRUE)

    expect_equal(dup2, dup)


})

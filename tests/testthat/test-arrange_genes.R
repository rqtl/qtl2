context("Arrange genes vertically so they don't overlap")

test_that("arrange_genes gives appropriate errors", {

    expect_error(arrange_genes(NULL, NULL))
    expect_error(arrange_genes(numeric(0), numeric(0)))
    expect_error(arrange_genes(c(1,2), 4))
    expect_error(arrange_genes(c(1,2), c(4,5,6)))

})

test_that("arrange_genes gives appropriate errors", {

    # no overlaps
    expect_equal(arrange_genes(1,2), 1)
    expect_equal(arrange_genes(c(1,5), c(2,6)), c(1,1))
    expect_equal(arrange_genes(c(1,5,10), c(2,6,11)), c(1,1,1))

    # some overlap
    expect_equal(arrange_genes(c(1,5), c(6,10)), c(1,2))
    expect_equal(arrange_genes(c(1,5,11), c(6,10,12)), c(1,2,1))
    expect_equal(arrange_genes(c(1,5,9), c(6,10,12)), c(1,2,1))
    expect_equal(arrange_genes(c(1,5,5.5), c(6,10,12)), c(1,2,3))

})

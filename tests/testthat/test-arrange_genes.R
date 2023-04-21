context("Arrange genes vertically so they don't overlap")

test_that("arrange_genes gives appropriate errors", {

#    expect_error(arrange_genes(NULL, NULL))
    expect_error(arrange_genes(numeric(0), numeric(0)))
    expect_error(arrange_genes(c(1,2), 4))
    expect_error(arrange_genes(c(1,2), c(4,5,6)))

})

test_that("arrange_genes works in simple cases", {

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


test_that("the example is arranged properly", {
    # this was a bug in my initial implementation
    # was giving 1, 1, 1, 1, 2, 1, 1, 2
    expect_equal( arrange_genes(c(139.99, 140.68, 141.71, 142.23, 142.59, 143.23, 144.40, 144.99),
                                c(140.77, 141.89, 142.74, 143.25, 143.48, 144.53, 145.46, 145.82)),

                  c(1, 2, 1, 2, 3, 1, 2, 1) )
})

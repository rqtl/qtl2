context("assign allele codes")

test_that("assign_allele_codes works", {

    expect_equal(assign_allele_codes(2, c("AA", "AB", "BB")), c("A", "B"))
    expect_equal(assign_allele_codes(2, c("AA", "AC", "CC")), c("A", "C"))

    expect_equal(assign_allele_codes(2), c("A", "B"))

    expect_equal(assign_allele_codes(3, c("AA", "AB", "BA", "BB", "AC", "CA", "CC", "BC", "CB")), c("A", "B", "C"))

})

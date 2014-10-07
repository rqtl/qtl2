
context("check_crosstype")

test_that("check_crosstype works appropriately", {

    good_types <- c("bc", "f2", "riself", "risib", "dh", "haploid")
    for(type in good_types)
        expect_true(check_crosstype(type))

    expect_error(check_crosstype("f2pk"))
    expect_false(check_crosstype("f2pk", FALSE))

    expect_error(check_crosstype("bada1gakadsf"))
    expect_false(check_crosstype("bada1gakadsf", FALSE))

})

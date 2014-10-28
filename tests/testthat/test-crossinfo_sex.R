context("check cross_info, sex, and handling of X")

test_that("checks of cross_info, sex, and X are correct", {

    null_crossinfo <- matrix(0L, 0, 0)
    null_isfemale <- logical(0)

    # dh, haploid, riself don't want X chr
    for(crosstype in c("dh", "haploid", "riself")) {
        expect_true(check_handle_x_chr(crosstype, FALSE))
        expect_false(check_handle_x_chr(crosstype, TRUE))
    }

    # others handle X chr ok
    for(crosstype in c("bc", "f2", "f2pk", "risib")) {
        expect_true(check_handle_x_chr(crosstype, FALSE))
        expect_true(check_handle_x_chr(crosstype, TRUE))
    }

    # dh, haploid, riself, risib ignore is_female
    for(crosstype in c("dh", "haploid", "riself", "risib")) {
        expect_true(check_is_female_vector(crosstype, null_isfemale, FALSE))
        expect_true(check_is_female_vector(crosstype, null_isfemale, TRUE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, FALSE), FALSE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, FALSE), TRUE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, NA, FALSE), FALSE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, NA, FALSE), TRUE))
    }

    # backcross and intercross need is_female
    for(crosstype in c("bc", "f2", "f2pk")) {
        expect_true(check_is_female_vector(crosstype, null_isfemale, FALSE))
        expect_false(check_is_female_vector(crosstype, null_isfemale, TRUE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, FALSE), TRUE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, FALSE), FALSE))
        expect_false(check_is_female_vector(crosstype, c(TRUE, NA, FALSE), TRUE))
        expect_true(check_is_female_vector(crosstype, c(TRUE, NA, FALSE), FALSE))
    }

    # dh, haploid, riself, bc ignore cross_info
    for(crosstype in c("dh", "haploid", "riself", "bc")) {
        expect_true(check_crossinfo(crosstype, null_crossinfo, FALSE))
        expect_true(check_crossinfo(crosstype, null_crossinfo, TRUE))
        expect_true(check_crossinfo(crosstype, cbind(c(0,1)), TRUE))
    }

    # intercross, phase-known intercross, and risib need cross_info for X
    for(crosstype in c("f2", "f2pk", "risib")) {
        expect_true(check_crossinfo(crosstype, null_crossinfo, FALSE))
        expect_false(check_crossinfo(crosstype, null_crossinfo, TRUE))
        expect_true(check_crossinfo(crosstype, cbind(c(0,1)), TRUE))
        expect_true(check_crossinfo(crosstype, cbind(c(0,1)), FALSE))
        expect_false(check_crossinfo(crosstype, cbind(c(0,NA,1)), TRUE))
        expect_true(check_crossinfo(crosstype, cbind(c(0,NA,1)), FALSE))
        expect_false(check_crossinfo(crosstype, cbind(c(0,2,1)), TRUE))
        expect_true(check_crossinfo(crosstype, cbind(c(0,2,1)), FALSE))
        expect_false(check_crossinfo(crosstype, cbind(0,0,1), TRUE))
        expect_true(check_crossinfo(crosstype, cbind(0,1,1), FALSE))
    }

})

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
    for(crosstype in c("bc", "f2", "f2pk", "risib", "ail")) {
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

    # backcross, intercross, and AIL need is_female
    for(crosstype in c("bc", "f2", "f2pk", "ail")) {
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

    # AIL needs no. generations + (for X chromosome) direction (0=AxB, 1=BxA, 2=balanced)
    expect_false(check_crossinfo("ail", null_crossinfo, FALSE))
    expect_false(check_crossinfo("ail", null_crossinfo, TRUE))
    expect_true(check_crossinfo("ail", cbind(c(2,3,25,50)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50)), TRUE))
    expect_true(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,2)), FALSE))
    expect_true(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,2)), TRUE))

    # no. generations >= 2 and not missing
    expect_false(check_crossinfo("ail", cbind(c(1,3,25,50),c(0,0,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,NA,25,50),c(0,0,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,1,50),c(0,0,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,NA),c(0,0,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(1,3,25,50),c(0,0,1,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,NA,25,50),c(0,0,1,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,1,50),c(0,0,1,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,NA),c(0,0,1,2)), TRUE))

    # dir = 0,1,2 and not missing
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,-1,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,-1)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(NA,0,1,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,NA,2)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,NA)), FALSE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,-1,1,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,-1)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(NA,0,1,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,NA,2)), TRUE))
    expect_false(check_crossinfo("ail", cbind(c(2,3,25,50),c(0,0,1,NA)), TRUE))

})

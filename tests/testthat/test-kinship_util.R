context("kinship utilities")

test_that("is_kinship and is_kinship_decomposed work", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[1:100,c(11,12,16,"X")]
    pr <- calc_genoprob(iron)

    k <- calc_kinship(pr)
    kloco <- calc_kinship(pr, "loco")

    kd <- decomp_kinship(k)
    kdloco <- decomp_kinship(kloco)

    expect_false(is_kinship_list(k))
    expect_false(is_kinship_list(kd))
    expect_true(is_kinship_list(kloco))
    expect_true(is_kinship_list(kdloco))

    expect_false(is_kinship_decomposed(k))
    expect_false(is_kinship_decomposed(kloco))
    expect_true(is_kinship_decomposed(kd))
    expect_true(is_kinship_decomposed(kdloco))

    expect_error(check_kinship_onechr(kloco))
    expect_error(check_kinship_onechr(kdloco))

    expect_equal(check_kinship_onechr(kloco[2]), kloco[[2]])
    expect_equivalent(check_kinship_onechr(kdloco[2]), kdloco[[2]]) # former will have the eigen_decomp attribute
    expect_equal(check_kinship_onechr(k), k)
    expect_equal(check_kinship_onechr(kd), kd)

})

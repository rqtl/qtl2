context("subset kinship matrices")

test_that("subset.kinship works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[c("2","3","4","5","6","7"), 2:6]
    pr <- calc_genoprob(iron)

    k <- calc_kinship(pr)
    k_loco <- calc_kinship(pr, "loco")

    # subset by individuals
    expect_equal(subset_kinship(k, ind=c(2,4)), k[c(2,4), c(2,4)])
    expect_equal(subset_kinship(k, ind=c("2","4")), k[c("2","4"), c("2","4")])

    # subset by individuals, loco case
    expect_equal(subset_kinship(k_loco, ind=c(2,4)), lapply(k_loco, function(a) a[c(2,4), c(2,4)]))
    expect_equal(subset_kinship(k_loco, ind=c("2","4")), lapply(k_loco, function(a) a[c("2","4"), c("2","4")]))

    # chr argument ignored if a plain matrix
    expect_equal(subset_kinship(k, ind=c(2,4), chr=5), k[c(2,4), c(2,4)])
    expect_equal(subset_kinship(k, ind=c("2","4"), chr=5), k[c("2","4"), c("2","4")])

    # chr argument, loco case, one chr in result
    expect_equal(subset_kinship(k_loco, ind=c(2,4), chr=5), k_loco[["5"]][c(2,4), c(2,4)])
    expect_equal(subset_kinship(k_loco, ind=c("2","4"), chr=5), k_loco[["5"]][c("2","4"), c("2","4")])

    # chr argument, loco case, mult chr in result
    expect_equal(subset_kinship(k_loco, ind=c(2,4), chr=c(2,5)), lapply(k_loco[c("2","5")], function(a) a[c(2,4), c(2,4)]))
    expect_equal(subset_kinship(k_loco, ind=c("2","4"), chr=c(2,5)), lapply(k_loco[c("2", "5")], function(a) a[c("2","4"), c("2","4")]))

    # chr argument, loco case, negative subscripts
    expect_equal(subset_kinship(k_loco, ind=c(2,4), chr=-c(2,4)), lapply(k_loco[c("3","5","6")], function(a) a[c(2,4), c(2,4)]))
    expect_equal(subset_kinship(k_loco, ind=c("2","4"), chr=c("-2","-4")), lapply(k_loco[c("3","5","6")], function(a) a[c("2","4"), c("2","4")]))

})

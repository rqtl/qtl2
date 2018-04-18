context("pull genotype probabilities for an interval")

test_that("pull_genoprobint works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(2,8,"X")]
    gmap <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, gmap, error_prob=0.002)

    pr_sub <- pull_genoprobint(pr, gmap, "8", c(25, 30))

    expected <- pr[,"8"]
    gm8 <- gmap[["8"]]
    markers <- names(gm8)[gm8 >= 25 & gm8 <= 30]
    expected[[1]] <- expected[[1]][,,markers,drop=FALSE]

    expect_equal(pr_sub, expected)

})

context("plot_genoprob")

test_that("plot_genoprob works", {

    skip_if(isnt_karl(), "plot tests only done locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[c("116", "232"),2]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    test_plot_genoprob <- function() {
        par(mfrow=c(2,1))
        plot_genoprob(probs, map, ind="116", main="116")
        plot_genoprob(probs, map, ind="232", main="232")
    }

    expect_doppelganger("plot_genoprob", test_plot_genoprob)

})

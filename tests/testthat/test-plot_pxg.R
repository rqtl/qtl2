context("plot_pxg")

test_that("plot_pxg works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    set.seed(62610474)
    geno <- maxmarg(probs, map, chr=16, pos=28.6, return_char=TRUE)

    test_plot_pxg <- function() plot_pxg(geno, log10(iron$pheno[,1]), ylab="log liver")
    vdiffr::expect_doppelganger("plot_pxg", test_plot_pxg)

    test_plot_pxg_se <- function() plot_pxg(geno, log10(iron$pheno[,1]), ylab="log liver", SEmult=2)
    vdiffr::expect_doppelganger("plot_pxg_se", test_plot_pxg_se)

})

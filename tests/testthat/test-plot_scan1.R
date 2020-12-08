context("plot_scan1")

test_that("plot_scan1 works", {

    skip_if(isnt_karl(), "plot tests only done locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)
    ylim <- c(0, maxlod(out)*1.02) # need to strip class to get overall max LOD
    chr <- c(2,7,8,9,15,16)

    test_plot_scan1 <- function() {
        plot(out, map, chr=chr, ylim=ylim)
        plot(out, map, lodcolumn=2, chr=chr, col="violetred", add=TRUE)
        legend("topleft", lwd=2, col=c("darkslateblue", "violetred"), colnames(out),
               bg="gray90") }

    expect_doppelganger("plot_scan1", test_plot_scan1)


    # single chromosome
    test_plot_scan1_onechr <- function() plot(out, map, chr=2)
    expect_doppelganger("plot_scan1_onechr", test_plot_scan1_onechr)

})

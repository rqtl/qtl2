context("plot_peaks")

test_that("plot_peaks works", {

    skip_if(isnt_karl(), "plot tests only done locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # find peaks above lod=3.5 (and calculate 1.5-LOD support intervals)
    peaks <- find_peaks(out, map, threshold=3.5, drop=1.5)

    test_plot_peaks <- function() plot_peaks(peaks, map)

    expect_doppelganger("plot_peaks", test_plot_peaks)

    peaks <- find_peaks(out, map, threshold=3.5)
    test_plot_peaks_noci <- function() plot_peaks(peaks, map)
    expect_doppelganger("plot_peaks_noci", test_plot_peaks_noci)

})

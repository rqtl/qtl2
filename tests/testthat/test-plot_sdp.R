context("plot_sdp")

test_that("plot_sdp works", {

    skip_if(isnt_karl(), "plot tests only run locally")

    pos <- c( 1,  6,  8,  9, 18, 22, 29, 38, 40, 48,
             57, 62, 69, 72, 74, 82, 83, 85, 90, 96)

    sdp <- c(127, 32, 177, 242, 248, 103, 188,  70,  76, 100,
             9,   19,  41,  13,  38, 202, 211, 105,  48, 135)

    test_plot_sdp <- function() plot_sdp(pos, sdp, xlim=c(0, 100))

    expect_doppelganger("plot_sdp", test_plot_sdp)

})

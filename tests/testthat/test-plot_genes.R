context("plot_genes")

test_that("plot_genes works", {

    skip_if(isnt_karl(), "plot tests only run locally")

    genes <- data.frame(chr = c("6", "6", "6", "6", "6", "6", "6", "6"),
                        start = c(139988753, 140680185, 141708118, 142234227, 142587862,
                                  143232344, 144398099, 144993835),
                        stop  = c(140041457, 140826797, 141773810, 142322981, 142702315,
                                  143260627, 144399821, 145076184),
                        strand = c("-", "+", "-", "-", "-", NA, "+", "-"),
                        Name = c("Plcz1", "Gm30215", "Gm5724", "Slco1a5", "Abcc9",
                                 "4930407I02Rik", "Gm31777", "Bcat1"),
                        stringsAsFactors=FALSE)

    test_plot_genes <- function() plot_genes(genes, xlim=c(140, 146), scale=1e-6)

    expect_doppelganger("plot_genes", test_plot_genes)

})

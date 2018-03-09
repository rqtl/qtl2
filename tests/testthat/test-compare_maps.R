context("compare maps")

test_that("compare_maps works", {

    iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
    gmap <- iron$gmap
    pmap <- iron$pmap

    expected <- data.frame(marker=character(0),
                           chr_map1=character(0),
                           pos_map1=numeric(0),
                           chr_map2=character(0),
                           pos_map2=numeric(0),
                           stringsAsFactors=FALSE)

    expect_equal(compare_maps(gmap, pmap), expected)

    # omit a marker from each map
    gmap[[7]] <- gmap[[7]][-3]
    pmap[[8]] <- pmap[[8]][-7]
    # swap order of a couple of markers on the physical map
    names(pmap[[9]])[3:4] <- names(pmap[[9]])[4:3]
    # move a marker to a different chromosome
    pmap[[10]] <- c(pmap[[10]], pmap[[1]][2])[c(1,3,2)]
    pmap[[1]] <- pmap[[1]][-2]

    expected <- structure(list(marker = c("D8Mit120", "D7Nds5", "D1Mit80", "D9Mit10", "D9Mit182"),
                               chr_map1 = c("8", NA, "1", "9", "9"),
                               pos_map1 = c(69.9, NA, 51.4, 43.7, 53.6),
                               chr_map2 = c(NA, "7", "10", "9", "9"),
                               pos_map2 = c(NA, 44.201096, 86.377953, 101.52887, 89.977078)),
                          .Names = c("marker", "chr_map1", "pos_map1", "chr_map2", "pos_map2"),
                          row.names = c(NA, 5L), class = "data.frame")

    expect_equal(compare_maps(gmap, pmap), expected)

})

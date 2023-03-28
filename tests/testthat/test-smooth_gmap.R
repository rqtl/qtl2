context("smooth_gmap")

test_that("smooth_gmap works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    # if alpha==0, nothing happens
    expect_equal(smooth_gmap(iron$gmap, iron$pmap, 0), iron$gmap)

    gmap_adj <- smooth_gmap(iron$gmap, iron$pmap, 0.02)

    # structure shouldn't change
    expect_equal(names(gmap_adj), names(iron$gmap))
    expect_equal(sapply(gmap_adj, length), sapply(iron$gmap, length))
    expect_equal(names(unlist(gmap_adj)), names(unlist(iron$gmap)))

    # chromosomes with just two markers are unchanged
    nmar <- sapply(gmap_adj, length)
    expect_equal( gmap_adj[nmar==2], iron$gmap[nmar==2] )

    # test unsmooth gmap
    gmap_back <- unsmooth_gmap(gmap_adj, iron$pmap, 0.02)
    expect_equal(iron$gmap, gmap_back)

    # test with a different alpha value
    expect_equal(unsmooth_gmap(smooth_gmap(iron$gmap, iron$pmap, 0.05), iron$pmap, 0.05), iron$gmap)

})

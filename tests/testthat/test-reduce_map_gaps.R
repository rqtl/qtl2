context("reduce gaps in a map")

test_that("reduce_map_gaps works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- iron$gmap

    expect_equal(reduce_map_gaps(map, 1000), map)

    expected <- map
    expected[[1]][3] <- expected[[1]][2] + 40
    expected[[4]][2] <- expected[[4]][1] + 40
    expected[[5]][2] <- expected[[5]][1] + 40

    expect_equal(reduce_map_gaps(map, 40), expected)

    map2 <- map
    map2[[2]][3:5] <- map2[[2]][3:5]+50
    expected[[2]][3:5] <- expected[[2]][3:5] - diff(expected[[2]][2:3]) + 40
    expect_equal(reduce_map_gaps(map2, 40), expected)

})

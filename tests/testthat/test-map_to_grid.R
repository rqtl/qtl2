context("Reduce map to pseudomarker grid")

test_that("map_to_grid works", {

    # try it out
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    grid <- calc_grid(grav2$gmap, step=1)

    reduced_map <- map
    for(i in seq(along=map)) {
        reduced_map[[i]] <- map[[i]][grid[[i]]]
    }

    expect_equal(reduced_map, map_to_grid(map, grid))
})

context("calculation of grid positions")
suppressMessages(library(qtl))

test_that("grid-based version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)
    names(map) <- paste0("map", seq(along=map))

    # step = marker distance
    pmap <- insert_pseudomarkers(list("1"=map), step=2.5, off_end=0, stepwidth="fixed")
    grid <- calc_grid(list("1"=map), step=2.5, off_end=0)
    expect_equal(length(pmap[[1]]), length(grid[[1]]))
    expect_equal(names(pmap[[1]]), names(grid[[1]]))

    expected_grid <- rep(TRUE, length(map))
    names(expected_grid) <- names(map)
    expect_equal(grid, list("1"=expected_grid))

    # step = 1
    pmap <- insert_pseudomarkers(list("1"=map), step=1, off_end=0, stepwidth="fixed")
    grid <- calc_grid(list("1"=map), step=1, off_end=0)

    # expected grid
    expected_grid <- !is.na(match(pmap[[1]], seq(0, 50, by=1)))
    names(expected_grid) <- names(pmap[[1]])
    expect_equal(grid, list("1"=expected_grid))

})

test_that("minimal version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)
    names(map) <- paste0("map", seq(along=map))

    # step = marker distance
    grid <- calc_grid(list("1"=map), step=2.5, off_end=0)
    expected_grid <- rep(TRUE, length(map))
    names(expected_grid) <- names(map)
    expect_equal(grid, list("1"=expected_grid))

})



test_that("insert_pseudomarkers gives distinct pseudomarker names with iron data", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    pmap <- insert_pseudomarkers(iron$gmap, step=1)
    grid <- calc_grid(iron$gmap, step=1)

    expect_equal(length(pmap), length(grid))
    expect_equal(names(pmap), names(grid))

    for(i in seq(along=pmap)) {
        expect_equal(length(pmap[[i]]), length(grid[[i]]))
        expect_equal(names(pmap[[i]]), names(grid[[i]]))
    }

})

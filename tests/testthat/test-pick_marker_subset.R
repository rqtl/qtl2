context("pick subset of markers")

test_that("pick_marker_subset matches qtl::pickMarkerSubset", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    map <- grav2$gmap

    wts <- lapply(vapply(map, length, 1),
                  function(n) runif(n, 1, 5))

    submap <- pick_marker_subset(map, 1, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- pick_marker_subset(map[[i]], 1, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 1, wts[[i]]))
    }

    # repeat with 5 cM
    submap <- pick_marker_subset(map, 5, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- pick_marker_subset(map[[i]], 5, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 5, wts[[i]]))
    }

    # without weights, get same result
    set.seed(8764972)
    submap <- pick_marker_subset(map, 1)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- pick_marker_subset(map[[i]], 1)
        expect_equal(submapi, submap[[i]])
    }

    # repeat at 5 cM
    set.seed(8764972)
    submap <- pick_marker_subset(map, 5)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- pick_marker_subset(map[[i]], 5)
        expect_equal(submapi, submap[[i]])
    }

})

context("create marker index")
suppressMessages(library(qtl))

test_that("grid-based version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)
    names(map) <- paste0("map", seq(along=map))

    # step = marker distance
    pmap <- insert_pseudomarkers(list("1"=map), step=2.5, off_end=0, stepwidth="fixed")
    index <- create_marker_index(list("1"=names(map)), pmap)
    expect_equal(index, list("1"=seq(along=map)-1))

    # step = 1
    pmap <- insert_pseudomarkers(list("1"=map), step=1, off_end=0, stepwidth="fixed")

    # expected index
    index <- rep(-1, length(pmap[[1]]))
    names(index) <- names(pmap[[1]])
    index[names(map)] <- seq(along=map)-1
    names(index) <- NULL

    expect_equal(list("1"=index), create_marker_index(list("1"=names(map)), pmap))
})

test_that("minimal version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)
    names(map) <- paste0("map", seq(along=map))

    # step = marker distance
    pmap <- insert_pseudomarkers(list("1"=map), step=2.5, off_end=0, stepwidth="max")
    expect_equal(list("1"=seq(along=map)-1), create_marker_index(list("1"=names(map)), pmap))

    # step = 1
    pmap <- insert_pseudomarkers(list("1"=map), step=1, off_end=0, stepwidth="max")

    # expected index
    index <- rep(-1, length(pmap[[1]]))
    names(index) <- names(pmap[[1]])
    index[names(map)] <- seq(along=map)-1
    names(index) <- NULL

    expect_equal(list("1"=index), create_marker_index(list("1"=names(map)), pmap))
})

test_that("minimal version works in more realistic case", {

    data(hyper)
    map <- qtl::pull.map(hyper, chr=1)

    pmap <- insert_pseudomarkers(map, step=1.55, off_end=4, stepwidth="max")

    expected_index <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
                        0, 0, 3, 0, 4, 0, 5, 0, 0, 6, 0, 7, 8, 0, 0, 0, 9, 0, 0, 0, 10,
                        0, 0, 0, 0, 0, 0, 11, 0, 0, 12, 0, 13, 0, 0, 14, 15, 0, 0, 0,
                        0, 16, 17, 18, 19, 0, 0, 20, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 22, 0, 0)-1
    expect_equal(list("1"=expected_index), create_marker_index(list("1"=names(map[[1]])), pmap))

})

test_that("insert_pseudomarkers works with a custom pseudomarker map", {

    data(hyper)
    map <- qtl::pull.map(hyper)

    set.seed(99735998)
    pseudomarker_map <- vector("list", length(map))
    names(pseudomarker_map) <- names(map)
    for(i in seq(along=map)) {
        n.pmar <- 10
        pseudomarker_map[[i]] <- sort(runif(n.pmar, 0, max(map[[i]])))
        names(pseudomarker_map[[i]]) <- paste0("c", names(map)[i], ".loc", 1:n.pmar)
    }

    combined_map <- insert_pseudomarkers(map, pseudomarker_map = pseudomarker_map)

    for(i in seq(along=map)) {
        expected <- match(names(combined_map[[i]]), names(map[[i]]))
        expected[is.na(expected)] <- 0
        expected <- expected - 1
        expect_equal(expected, create_marker_index(list("1"=names(map[[i]])), list("1"=combined_map[[i]]))[[1]])
    }
})

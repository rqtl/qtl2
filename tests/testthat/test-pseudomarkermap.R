
context("construction of pseudomarker map")
library(qtl)

test_that("grid-based version works in simple case", {

    #equally-spaced map
    map <- seq(0, 50, by=2.5)

    # step = marker distance
    pmap <- create_pseudomarker_map_grid(map, step=2.5, off_end=0)
    expect_equivalent(map, pmap)
    expect_equal(attr(pmap, "index"), seq(along=map))
    expect_equal(attr(pmap, "grid"), rep(TRUE, length(map)))

    # step = 1
    pmap <- create_pseudomarker_map_grid(map, step=1, off_end=0)
    pmap_qtl <- qtl::create.map(map, step=1, off.end=0)
    expect_equivalent(pmap, pmap_qtl)

    # expected index
    index <- rep(0, length(pmap))
    mar <- names(pmap)==""
    index[mar] <- 1:sum(mar)
    expect_equal(attr(pmap, "index"), index)

    # expected grid
    expect_equal(attr(pmap, "grid"), !is.na(match(pmap, seq(0, 50, by=1))))

})

# the following are totally not working
# ...the new routine seems better than the old one
test_that("grid-based version works in more realistic case", {

    data(hyper)

    map <- qtl::pull.map(hyper, chr=1)[[1]]

    pmap <- create_pseudomarker_map_grid(map, step=1.55, off_end=4)

    expected <- c(-0.7, 0.85, 2.4, 3.3, 3.95, 5.5, 7.05, 8.6, 10.15,
                  11.7, 13.25, 14.8, 16.35, 17.9, 19.45, 19.7000000001, 21, 22.55,
                  24.1, 25.65, 27.2, 28.75, 30.3, 31.85, 32.8000000002, 33.4, 34.95,
                  35.0000000003, 36.5, 37.2000000004, 38.05, 39.6, 41.15, 41.5000000005,
                  42.7, 43.7000000006, 43.7000000007, 44.25, 45.8, 47.35, 48.9,
                  49.2000000008, 50.45, 52, 53.55, 54.6000000009, 55.1, 56.65,
                  58.2, 59.75, 61.3, 62.85, 64.4, 64.500000001, 65.95, 67.5, 67.8000000011,
                  69.05, 69.9000000012, 70.6, 72.15, 73.7, 74.3000000013, 75.25,
                  75.4000000014, 76.8, 78.35, 79.9, 81.45, 82.0000000015, 82.0000000016,
                  82.0000000017, 82.0000000018, 83, 84.55, 86.1, 86.3000000019,
                  87.65, 89.2, 90.75, 92.3, 93.85, 94.000000002, 95.4, 96.95, 98.5,
                  100.05, 101.6, 103.15, 104.7, 106.25, 107.8, 109.35, 110.9, 112.45,
                  114, 115.55, 115.8000000021, 117.1, 118.65)
    pmap_names <- c("loc-1", "loc1", "loc2", "D1Mit296", "loc4", "loc6", "loc7", "loc9", "loc10",
                    "loc12", "loc13", "loc15", "loc16", "loc18", "loc19", "D1Mit123",
                    "loc21", "loc23", "loc24", "loc26", "loc27", "loc29", "loc30",
                    "loc32", "D1Mit156", "loc33", "loc35", "D1Mit178", "loc36", "D1Mit19",
                    "loc38", "loc40", "loc41", "D1Mit7", "loc43", "D1Mit46", "D1Mit132",
                    "loc44", "loc46", "loc47", "loc49", "D1Mit334", "loc50", "loc52",
                    "loc54", "D1Mit305", "loc55", "loc57", "loc58", "loc60", "loc61",
                    "loc63", "loc64", "D1Mit26", "loc66", "loc68", "D1Mit94", "loc69",
                    "D1Mit218", "loc71", "loc72", "loc74", "D1Mit100", "loc75", "D1Mit102",
                    "loc77", "loc78", "loc80", "loc81", "D1Mit14", "D1Mit105", "D1Mit159",
                    "D1Mit267", "loc83", "loc85", "loc86", "D1Mit15", "loc88", "loc89",
                    "loc91", "loc92", "loc94", "D1Mit456", "loc95", "loc97", "loc98",
                    "loc100", "loc102", "loc103", "loc105", "loc106", "loc108", "loc109",
                    "loc111", "loc112", "loc114", "loc116", "D1Mit155", "loc117",
                    "loc119")

    expect_equivalent(pmap, expected)
    expect_equal(names(pmap), pmap_names)

    expected_index <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
                        0, 0, 0, 0, 3, 0, 0, 4, 0, 5, 0, 0, 0, 6, 0, 7, 8, 0, 0, 0, 0,
                        9, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 12, 0, 13, 0,
                        0, 0, 14, 0, 15, 0, 0, 0, 0, 16, 17, 18, 19, 0, 0, 0, 20, 0,
                        0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 22,
                        0, 0)
    expect_equal(attr(pmap, "index"), expected_index)

    expected_grid <- c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                       TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                       TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE,
                       TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE,
                       TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                       TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE,
                       TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE,
                       FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
                       TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
    expect_equal(attr(pmap, "grid"), expected_grid)


})


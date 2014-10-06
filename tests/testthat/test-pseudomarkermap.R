
context("construction of pseudomarker map")
library(qtl)

test_that("grid-based version works in simple case", {

    #equally-spaced map
    map <- seq(0, 50, by=2.5)

    # step = marker distance
    pmap <- insert_pseudomarkers(map, step=2.5, off_end=0, stepwidth="fixed")
    expect_equivalent(map, pmap)
    expect_equal(attr(pmap, "index"), seq(along=map))
    expect_equal(attr(pmap, "grid"), rep(TRUE, length(map)))

    # step = 1
    pmap <- insert_pseudomarkers(map, step=1, off_end=0, stepwidth="fixed")
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

test_that("minimal version works in simple case", {

    #equally-spaced map
    map <- seq(0, 50, by=2.5)

    # step = marker distance
    pmap <- insert_pseudomarkers(map, step=2.5, off_end=0, stepwidth="max")
    expect_equivalent(map, pmap)
    expect_equal(attr(pmap, "index"), seq(along=map))

    # step = 1
    pmap <- insert_pseudomarkers(map, step=1, off_end=0, stepwidth="max")
    expect_equivalent(pmap, seq(0, 50, by=5/6))

    # expected index
    index <- rep(0, length(pmap))
    mar <- names(pmap)==""
    index[mar] <- 1:sum(mar)
    expect_equal(attr(pmap, "index"), index)

})

# the following are totally not working
# ...the new routine seems better than the old one
test_that("minimal version works in more realistic case", {

    data(hyper)
    map <- qtl::pull.map(hyper, chr=1)[[1]]

    pmap <- insert_pseudomarkers(map, step=1.55, off_end=4, stepwidth="max")

    expected <- c(0.2, 1.75, 3.3, 4.79090909091818, 6.28181818183636, 7.77272727275454,
                  9.26363636367273, 10.7545454545909, 12.2454545455091, 13.7363636364273,
                  15.2272727273455, 16.7181818182636, 18.2090909091818, 19.7000000001,
                  21.1555555556667, 22.6111111112333, 24.0666666668, 25.5222222223667,
                  26.9777777779333, 28.4333333335, 29.8888888890667, 31.3444444446333,
                  32.8000000002, 33.90000000025, 35.0000000003, 36.10000000035,
                  37.2000000004, 38.6333333337667, 40.0666666671333, 41.5000000005,
                  42.60000000055, 43.7000000006, 43.7000000007, 45.075000000725,
                  46.45000000075, 47.825000000775, 49.2000000008, 50.550000000825,
                  51.90000000085, 53.250000000875, 54.6000000009, 56.0142857152,
                  57.4285714295, 58.8428571438, 60.2571428581, 61.6714285724, 63.0857142867,
                  64.500000001, 65.6000000010333, 66.7000000010667, 67.8000000011,
                  68.85000000115, 69.9000000012, 71.3666666679, 72.8333333346,
                  74.3000000013, 75.4000000014, 76.72000000142, 78.04000000144,
                  79.36000000146, 80.68000000148, 82.0000000015, 82.0000000016,
                  82.0000000017, 82.0000000018, 83.4333333351667, 84.8666666685333,
                  86.3000000019, 87.84000000192, 89.38000000194, 90.92000000196,
                  92.46000000198, 94.000000002, 95.45333333534, 96.90666666868,
                  98.36000000202, 99.81333333536, 101.2666666687, 102.72000000204,
                  104.17333333538, 105.62666666872, 107.08000000206, 108.5333333354,
                  109.98666666874, 111.44000000208, 112.89333333542, 114.34666666876,
                  115.8000000021, 117.3500000021, 118.9000000021)
    pmap_names <- c("loc0", "loc2", "D1Mit296", "loc5", "loc6", "loc8", "loc9",
                    "loc11", "loc12", "loc14", "loc15", "loc17", "loc18", "D1Mit123",
                    "loc21", "loc23", "loc24", "loc26", "loc27", "loc28", "loc30",
                    "loc31", "D1Mit156", "loc34", "D1Mit178", "loc36", "D1Mit19",
                    "loc39", "loc40", "D1Mit7", "loc43", "D1Mit46", "D1Mit132", "loc45",
                    "loc46", "loc48", "D1Mit334", "loc51", "loc52", "loc53", "D1Mit305",
                    "loc56", "loc57", "loc59", "loc60", "loc62", "loc63", "D1Mit26",
                    "loc66", "loc67", "D1Mit94", "loc69", "D1Mit218", "loc71", "loc73",
                    "D1Mit100", "D1Mit102", "loc77", "loc78", "loc79", "loc81", "D1Mit14",
                    "D1Mit105", "D1Mit159", "D1Mit267", "loc83", "loc85", "D1Mit15",
                    "loc88", "loc89", "loc91", "loc92", "D1Mit456", "loc95", "loc97",
                    "loc98", "loc100", "loc101", "loc103", "loc104", "loc106", "loc107",
                    "loc109", "loc110", "loc111", "loc113", "loc114", "D1Mit155",
                    "loc117", "loc119")

    expect_equivalent(pmap, expected)
    expect_equal(names(pmap), pmap_names)

    expected_index <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
                        0, 0, 3, 0, 4, 0, 5, 0, 0, 6, 0, 7, 8, 0, 0, 0, 9, 0, 0, 0, 10,
                        0, 0, 0, 0, 0, 0, 11, 0, 0, 12, 0, 13, 0, 0, 14, 15, 0, 0, 0,
                        0, 16, 17, 18, 19, 0, 0, 20, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 22, 0, 0)
    expect_equal(attr(pmap, "index"), expected_index)

})


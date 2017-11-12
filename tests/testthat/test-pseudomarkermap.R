context("construction of pseudomarker map")
suppressMessages(library(qtl))

test_that("insert_pseudomarkers gives an error if the input is NULL", {
    expect_error(insert_pseudomarkers(NULL, step=1))

    library(qtl)
    data(hyper)
    expect_error(insert_pseudomarkers(hyper$gmap, step=1, error.prob=0.01))

})


test_that("grid-based version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)

    # step = marker distance
    pmap <- insert_pseudomarkers(list("1"=map), step=2.5, off_end=0, stepwidth="fixed")
    expect_equivalent(map, pmap[[1]])

    # step = 1
    pmap <- insert_pseudomarkers(list("1"=map), step=1, off_end=0, stepwidth="fixed")
    pmap_qtl <- qtl::create.map(map, step=1, off.end=0)
    expect_equivalent(pmap[[1]], pmap_qtl)

})

test_that("minimal version works in simple case", {

    # equally-spaced map
    map <- seq(0, 50, by=2.5)

    # step = marker distance
    pmap <- insert_pseudomarkers(list("1"=map), step=2.5, off_end=0, stepwidth="max")
    expect_equivalent(map, pmap[[1]])

    # step = 1
    pmap <- insert_pseudomarkers(list("1"=map), step=1, off_end=0, stepwidth="max")
    expect_equivalent(pmap[[1]], seq(0, 50, by=5/6))

})

test_that("minimal version works in more realistic case", {

    data(hyper)
    map <- qtl::pull.map(hyper, chr=1)

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
    pmap_names <- c("c1.loc0", "c1.loc2", "D1Mit296", "c1.loc5", "c1.loc6", "c1.loc8", "c1.loc9",
                    "c1.loc11", "c1.loc12", "c1.loc14", "c1.loc15", "c1.loc17", "c1.loc18", "D1Mit123",
                    "c1.loc21", "c1.loc23", "c1.loc24", "c1.loc26", "c1.loc27", "c1.loc28", "c1.loc30",
                    "c1.loc31", "D1Mit156", "c1.loc34", "D1Mit178", "c1.loc36", "D1Mit19",
                    "c1.loc39", "c1.loc40", "D1Mit7", "c1.loc43", "D1Mit46", "D1Mit132", "c1.loc45",
                    "c1.loc46", "c1.loc48", "D1Mit334", "c1.loc51", "c1.loc52", "c1.loc53", "D1Mit305",
                    "c1.loc56", "c1.loc57", "c1.loc59", "c1.loc60", "c1.loc62", "c1.loc63", "D1Mit26",
                    "c1.loc66", "c1.loc67", "D1Mit94", "c1.loc69", "D1Mit218", "c1.loc71", "c1.loc73",
                    "D1Mit100", "D1Mit102", "c1.loc77", "c1.loc78", "c1.loc79", "c1.loc81", "D1Mit14",
                    "D1Mit105", "D1Mit159", "D1Mit267", "c1.loc83", "c1.loc85", "D1Mit15",
                    "c1.loc88", "c1.loc89", "c1.loc91", "c1.loc92", "D1Mit456", "c1.loc95", "c1.loc97",
                    "c1.loc98", "c1.loc100", "c1.loc101", "c1.loc103", "c1.loc104", "c1.loc106", "c1.loc107",
                    "c1.loc109", "c1.loc110", "c1.loc111", "c1.loc113", "c1.loc114", "D1Mit155",
                    "c1.loc117", "c1.loc119")

    expect_equivalent(pmap[[1]], expected)
    expect_equal(names(pmap[[1]]), pmap_names)
})

test_that("insert_pseudomarkers works with a custom pseudomarker map", {

    data(hyper)
    map <- qtl::pull.map(hyper)

    set.seed(99735998)
    pseudomarker_map <- vector("list", length(map))
    for(i in seq(along=map)) {
        n.pmar <- 10
        pseudomarker_map[[i]] <- sort(runif(n.pmar, 0, max(map[[i]])))
        names(pseudomarker_map[[i]]) <- paste0("c", names(map)[i], ".loc", 1:n.pmar)
    }

    combined_map <- insert_pseudomarkers(map, pseudomarker_map = pseudomarker_map)

    expect_equal(sapply(map, length) + sapply(pseudomarker_map, length),
                 sapply(combined_map, length))

    for(i in seq(along=map)) {
        combined_map_2 <- sort(c(map[[i]], pseudomarker_map[[i]]))
        expect_equivalent(combined_map_2, combined_map[[i]])
        expect_equal(names(combined_map_2), names(combined_map[[i]]))
    }

    expect_equal(names(map), names(combined_map))

})

test_that("insert_pseudomarkers gives distinct pseudomarker names with iron data", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    pmap <- insert_pseudomarkers(iron$gmap, step=1)

    expect_true(all(vapply(pmap$map, function(a) !any_duplicates(names(a)), TRUE)))

})

test_that("insert_pseudomarkers works with multi-core", {

    if(isnt_karl()) skip("this test only run locally")

    data(hyper)
    map <- qtl::pull.map(hyper)

    set.seed(99735998)
    pseudomarker_map <- vector("list", length(map))
    for(i in seq(along=map)) {
        n.pmar <- 10
        pseudomarker_map[[i]] <- sort(runif(n.pmar, 0, max(map[[i]])))
        names(pseudomarker_map[[i]]) <- paste0("c", names(map)[i], ".loc", 1:n.pmar)
    }

    combined_map <- insert_pseudomarkers(map, pseudomarker_map = pseudomarker_map, cores=8)

    expect_equal(sapply(map, length) + sapply(pseudomarker_map, length),
                 sapply(combined_map, length))

    for(i in seq(along=map)) {
        combined_map_2 <- sort(c(map[[i]], pseudomarker_map[[i]]))
        expect_equivalent(combined_map_2, combined_map[[i]])
        expect_equal(names(combined_map_2), names(combined_map[[i]]))
    }

    expect_equal(names(map), names(combined_map))

})

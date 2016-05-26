context("find_marker and find_markerpos")

test_that("find_marker works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # find markers by their genetic map positions
    expect_equal(find_marker(iron$gmap, c(8, 11), c(37.7, 56.9)), c("D8Mit294", "D11Mit101"))

    # find markers by their physical map positions (two markers on chr 7)
    expect_equal(find_marker(iron$pmap, 7, c(44.2, 108.9)), c("D7Nds5", "D7Mit17"))

    # find markers in an interval
    expect_equal(find_marker(iron$pmap, 16, interval=c(35, 80)), c("D16Mit4", "D16Mit30", "D16Mit19"))

    # error if you give both pos and interval
    expect_error(find_marker(iron$gmap, 17, 20, c(5,6)))

    # error if chr and pos are different lengths > 1
    expect_error(find_marker(iron$gmap, c(2,17), c(20, 40, 50)))

    # error if length(interval) != 2
    expect_error(find_marker(iron$pmap, 18, interval=20))
    expect_error(find_marker(iron$pmap, 18, interval=c(20, 40, 50)))

    # error if interval provided and length(chr) != 1
    expect_error(find_marker(iron$gmap, c(18, 19), interval=c(20, 30)))

})

test_that("find_markerpos works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    expected <- data.frame(chr=c("8", "11"),
                           gmap=c(39.1, 56.9),
                           pmap=c(82.491410, 112.653628),
                           stringsAsFactors=FALSE)
    rownames(expected) <- c("D8Mit294", "D11Mit101")
    expect_equal(find_markerpos(iron, c("D8Mit294", "D11Mit101")), expected)

    # not-found marker with na.rm=TRUE
    expect_equal(suppressWarnings(find_markerpos(iron, c("D8Mit294", "blah", "D11Mit101"))), expected)
    expect_warning(find_markerpos(iron, c("D8Mit294", "blah", "D11Mit101")))

    # not-found marker with na.rm=FALSE
    expected <- rbind(expected[1,], "blah"=data.frame(chr=NA, gmap=NA, pmap=NA), expected[2,])
    expect_equal(find_markerpos(iron, c("D8Mit294", "blah", "D11Mit101"), na.rm=FALSE), expected)

    # one marker
    expect_equal(find_markerpos(iron, "D8Mit294"), expected[1,])

    # no gmap
    iron_nogmap <- iron
    iron_nogmap$gmap <- NULL
    expect_equal(find_markerpos(iron_nogmap, c("D8Mit294", "D11Mit101")), expected[c(1,3),-2])
    expect_equal(find_markerpos(iron_nogmap, "D11Mit101"), expected[3,-2])

    # no pmap
    iron_nopmap <- iron
    iron_nopmap$pmap <- NULL
    expect_equal(find_markerpos(iron_nopmap, c("D8Mit294", "D11Mit101")), expected[c(1,3),-3])
    expect_equal(find_markerpos(iron_nopmap, "D11Mit101"), expected[3,-3])

    # no map at all
    iron_nomap <- iron_nopmap
    iron_nomap$gmap <- NULL
    expect_error(find_markerpos(iron_nomap, "D11Mit101"))

})

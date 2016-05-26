context("find_marker")

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

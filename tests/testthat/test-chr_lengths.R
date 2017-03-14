context("chromosome lengths")


test_that("chr_lengths works", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    expected <- c("1"=83.1, "2"=34.9, "3"=29.5, "4"=42.7, "5"=44.8,
                  "6"=25.2, "7"=52.5, "8"=75.3, "9"=54.6, "10"=33.9,
                  "11"=56.9, "12"=38.2, "13"=22.9, "14"=32.8, "15"=32.8,
                  "16"=44.8, "17"=36.0, "18"=26.2, "19"=35.0, "X"=28.4)
    is_x_chr <- (names(expected)=="X")
    names(is_x_chr) <- names(expected)
    attr(expected, "is_x_chr") <- is_x_chr

    expect_equal(chr_lengths(iron$gmap), expected)

    expected_collapsed <- c(A=802.1, X=28.4)
    attr(expected_collapsed, "is_x_chr") <- c(A=FALSE, X=TRUE)

    expect_equal(chr_lengths(iron$gmap, TRUE), expected_collapsed)

    expect_equal(collapse_chr_lengths_to_AX(expected), expected_collapsed)

    is_x_chr <- attr(expected, "is_x_chr")
    attr(expected, "is_x_chr") <- NULL

    expect_warning(collapse_chr_lengths_to_AX(expected))

    expect_equal(collapse_chr_lengths_to_AX(expected, is_x_chr), expected_collapsed)

})

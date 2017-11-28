context("chromosome lengths")


test_that("chr_lengths works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    # regular case
    expected <- c("1"=83.1, "2"=34.9, "3"=29.5, "4"=42.7, "5"=44.8,
                  "6"=25.2, "7"=52.5, "8"=75.3, "9"=54.6, "10"=33.9,
                  "11"=56.9, "12"=38.2, "13"=22.9, "14"=32.8, "15"=32.8,
                  "16"=44.8, "17"=36.0, "18"=26.2, "19"=35.0, "X"=28.4)
    is_x_chr <- (names(expected)=="X")
    names(is_x_chr) <- names(expected)
    attr(expected, "is_x_chr") <- is_x_chr
    expect_equal(chr_lengths(iron$gmap), expected)

    # check collapsing -> A, X lengths
    expected_collapsed <- c(A=802.1, X=28.4)
    attr(expected_collapsed, "is_x_chr") <- c(A=FALSE, X=TRUE)
    expect_equal(chr_lengths(iron$gmap, TRUE), expected_collapsed)
    expect_equal(collapse_chr_lengths_to_AX(expected), expected_collapsed)

    # problems with is_x_chr
    is_x_chr <- attr(expected, "is_x_chr")
    attr(expected, "is_x_chr") <- NULL
    expect_warning(collapse_chr_lengths_to_AX(expected))
    expect_equal(collapse_chr_lengths_to_AX(expected, is_x_chr), expected_collapsed)

    # already collapsed
    pre_collapsed <- c(A=1300, X=100)
    expected <- pre_collapsed
    attr(expected, "is_x_chr") <- c(A=FALSE, X=TRUE)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(FALSE,FALSE)), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(TRUE,TRUE)), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(TRUE,FALSE)), expected)

    # already collapsed, no names
    names(pre_collapsed) <- NULL
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(FALSE,FALSE)), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(TRUE,TRUE)), expected)
    expect_equal(collapse_chr_lengths_to_AX(pre_collapsed, c(TRUE,FALSE)), expected)

    # length two but not collapsed, really
    only_A <- pre_collapsed
    names(only_A) <- c("I", "II")
    attr(only_A, "is_x_chr") <- c(I=FALSE, II=FALSE)
    expected <- c("A"=1400, "X"=0)
    attr(expected, "is_x_chr") <- c(A=FALSE, X=TRUE)
    expect_equal(collapse_chr_lengths_to_AX(only_A), expected)
    expect_equal(collapse_chr_lengths_to_AX(only_A, c(FALSE,FALSE)), expected)
    expected[1] <- 0; expected[2] <- 1400
    expect_equal(collapse_chr_lengths_to_AX(only_A, c(TRUE,TRUE)), expected)
    expected[1] <- 100; expected[2] <- 1300
    expect_equal(collapse_chr_lengths_to_AX(only_A, c(TRUE,FALSE)), expected)

})

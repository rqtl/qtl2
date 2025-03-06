context("compare_founder_geno")

test_that("compare_founder_geno gives error with no founders", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    expect_error(compare_founder_geno(iron))

})


test_that("compare_founder_geno works", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    cg <- compare_founder_geno(DOex)

    expected <- structure(c(320, 320, 319, 319, 319, 308, 316, 318, 0.54375,
                            322, 321, 320, 321, 308, 316, 319, 0.595611285266458, 0.548286604361371,
                            321, 319, 321, 308, 315, 318, 0.67398119122257, 0.5625, 0.623824451410658,
                            320, 319, 307, 315, 317, 0.623824451410658, 0.526479750778816,
                            0.604361370716511, 0.61128526645768, 321, 308, 315, 318, 0.376623376623377,
                            0.301948051948052, 0.386363636363636, 0.374592833876221, 0.353896103896104,
                            308, 306, 306, 0.382911392405063, 0.268987341772152, 0.387301587301587,
                            0.349206349206349, 0.396825396825397, 0.679738562091503, 316,
                            314, 0.405660377358491, 0.347962382445141, 0.405660377358491,
                            0.384858044164038, 0.437106918238994, 0.594771241830065, 0.589171974522293,
                            319), dim = c(8L, 8L), dimnames = list(c("A", "B", "C", "D", "E", "F", "G", "H"),
                                                                   c("A", "B", "C", "D", "E", "F", "G", "H")),
                          proportion = TRUE, class = c("compare_geno", "matrix"))

    expect_equal(cg, expected)

})

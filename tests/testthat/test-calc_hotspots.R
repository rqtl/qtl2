context("calc hotspots")

test_that("calc_hotspots works", {

    skip_on_cran()

    url <- "https://kbroman.org/qtl2/assets/sampledata/pqtl_data.RData"
    tempfile <- file.path(tempdir(), basename(url))
    download.file(url, tempfile)
    load(tempfile)
    unlink(tempfile)

    map <- map[c(7,11)]
    map[["7"]] <- map[["7"]][round(seq(1, length(map[["7"]]), length=21))]
    map[["11"]] <- map[["11"]][round(seq(1, length(map[["11"]]), length=21))]

    result <- calc_hotspots(qtl, map, window=2)

    expected <- structure(c(2L, 3L, 1L, 5L, 2L, 2L, 0L, 1L, 3L, 2L, 0L, 7L, 89L,
                            9L, 1L, 3L, 3L, 3L, 5L, 4L, 2L, 2L, 2L, 4L, 1L, 33L, 3L, 8L,
                            5L, 10L, 5L, 1L, 4L, 2L, 1L, 10L, 1L, 3L, 4L, 9L, 7L, 1L),
                          dim = c(42L, 1L),
                          dimnames = list(c("ICR1775", "UNCHS019493", "UNCHS019598",
                                            "UNC12613387", "UNC12680825", "UNCHS019924", "JAX00640534", "JAX00152385",
                                            "UNC13107446", "UNC13214549", "UNCHS020578", "UNC13432293", "ICR5875",
                                            "UNC13623993", "JAX00652501", "UNCHS021418", "UNC13809062", "UNC13876294",
                                            "JAX00157485", "JAX00658243", "UNC14050946", "UNCJPD004490",
                                            "UNC19047951", "UNCHS029895", "UNCJPD004551", "UNCHS030186",
                                            "UNCHS030349", "JAX00309706", "UNCHS030603", "UNCHS030775", "UNC19754812",
                                            "UNCHS031128", "UNC19898072", "UNCHS031478", "UNCHS031633", "JAX00317756",
                                            "UNCHS031954", "JAX00320253", "UNC20324233", "UNC20412654", "UNC20479736",
                                            "UNCJPD004845"), "num_qtl"),
                          class = c("scan1", "matrix"))

    expect_equal(result, expected)
})

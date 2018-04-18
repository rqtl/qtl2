context("find index snp")

test_that("find_index_snp works", {

    snpinfo <- structure(list(chr = c("19", "19", "X", "X"),
                              pos = c(40.36, 40.53, 110.91, 111.21), sdp = c(170, 170, 240, 170),
                              snp = c("m1", "m2", "m3", "m4"), index = c(1L, 1L, 1L, 2L),
                              interval = c(91L, 91L, 236L, 236L),
                              on_map = c(FALSE, FALSE, FALSE, FALSE)),
                         .Names = c("chr", "pos", "sdp", "snp", "index", "interval", "on_map"),
                         row.names = c("m1", "m2", "m3", "m4"), class = "data.frame")

    expect_equal(find_index_snp(snpinfo, "m2"), "m1")
    expect_warning(expect_equal(find_index_snp(snpinfo, "blah"), as.character(NA)))
    expect_equal(find_index_snp(snpinfo, "m1"), "m1")
    expect_equal(find_index_snp(snpinfo, c("m1", "m2", "m3", "m4")), c("m1", "m1", "m3", "m4"))

})

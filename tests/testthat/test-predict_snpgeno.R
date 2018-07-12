context("predict_snpgeno")

test_that("predict_snpgeno works", {

    if(isnt_karl()) skip("this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)
    m <- maxmarg(probs)

    infg <- predict_snpgeno(DOex, m)
    expect_equal(class(infg), "list")
    expect_equal(length(infg), 1)
    expect_equal(names(infg), "2")
    expect_equal(dim(infg[[1]]), c(20,125))
    expect_true(all(is.na(infg[[1]][,1])))
    expect_equivalent(infg[[1]][9,1:10], c(NA, 2, 2, 1, 1, 2, 2, NA, NA, NA))
    expect_equivalent(infg[[1]][10,101:110], c(2,1,2,3,1,2,1,1,2,2))
    expect_equivalent(infg[[1]][11,51:60], c(NA,NA,NA,NA,NA,NA,3,2,2,NA))

})

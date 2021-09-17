context("predict_snpgeno")

test_that("predict_snpgeno works", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)
    m <- maxmarg(probs)

    infg <- predict_snpgeno(DOex, m)
    expect_equal(class(infg), "list")
    expect_equal(length(infg), 1)
    expect_equal(names(infg), "2")
    expect_equal(dim(infg[[1]]), c(20,127))
    expect_true(all(is.na(infg[[1]][,1])))
    expect_equivalent(infg[[1]][9,1:10], c(NA, 2, 1, 1, 3, 2, 2, NA, NA, NA))
    expect_equivalent(infg[[1]][10,101:110], c(2,2,2,2,2,2,2,1,3,2))
    expect_equivalent(infg[[1]][11,51:60], c(NA,NA,NA,NA,1,3,1,2,2,3))

})


test_that("predict_snpgeno works for magic lines", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/ArabMAGIC/arabmagic_tair9.zip")
    magic <- read_cross2(file)
    ind <- paste0("MAGIC", ".", 1:20)
    chr <- "2"
    probs <- calc_genoprob(magic[ind,chr], error_prob=0.002)
    m <- maxmarg(probs)

    infg <- predict_snpgeno(magic, m)

    expect_equal(class(infg), "list")
    expect_equal(length(infg), 1)
    expect_equal(names(infg), chr)
    expect_equal(dim(infg[[1]]), c(20,211))
    expect_equivalent(infg[[1]][,1], c(3,1,3,1,NA,3,1,3,1,1,3,3,NA,1,3,NA,1,1,NA,1))
    expect_equivalent(infg[[1]][9,1:10], c(1,1,1,3,3,1,1,3,3,1))
    expect_equivalent(infg[[1]][10,101:110],  c(rep(1,9), 3))
    expect_equivalent(infg[[1]][11,51:60], c(1,1,3,1,1,1,1,1,1,1))

    fg <- magic$founder_geno[[chr]]
    fg[fg==0] <- NA
    infg_hard <- t(sapply(1:20, function(wh_ind) sapply(seq_along(m[[chr]][wh_ind,]), function(i) {
        g <- m[[chr]][wh_ind,i]; ifelse(is.na(g), NA, fg[g,i]) })))
    expect_equivalent(infg[[chr]], infg_hard)

})

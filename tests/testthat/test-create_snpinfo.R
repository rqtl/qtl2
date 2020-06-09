context("create snpinfo from cross2 object")

test_that("create_snpinfo works", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")

    DOex <- read_cross2(file)

    snpinfo <- create_snpinfo(DOex)

    fg <- do.call("cbind", DOex$founder_geno)

    expect_equal(sum(colSums(fg==0)==0), nrow(snpinfo))

    n <- sapply(DOex$founder_geno, function(a) sum(colSums(a==0)==0))
    expect_equivalent(unclass(table(snpinfo$chr)), n)

    expect_equivalent(snpinfo$sdp, calc_sdp(t(fg[,colSums(fg==0)==0])))

})

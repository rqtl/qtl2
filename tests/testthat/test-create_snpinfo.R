context("create snpinfo from cross2 object")

test_that("create_snpinfo works", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")

    DOex <- read_cross2(file)

    snpinfo <- create_snpinfo(DOex)

    fg <- do.call("cbind", DOex$founder_geno)

    expect_equal(sum(colSums(fg==0)==0), nrow(snpinfo))

    n <- sapply(DOex$founder_geno, function(a) sum(colSums(a==0)==0))
    expect_equivalent(unclass(table(snpinfo$chr)), n)

    expect_equivalent(snpinfo$sdp, calc_sdp(t(fg[,colSums(fg==0)==0])))

})


test_that("create_snpinfo works for arab data", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/main/ArabMAGIC/arabmagic_tair9.zip")
    arab <- read_cross2(file)

    snpinfo <- create_snpinfo(arab)

    # expected markers, dropping non-informative or with missing founder genotype
    tmp <- drop_nullmarkers(arab)
    fg <- do.call("cbind", tmp$founder_geno)
    fg <- fg[,colSums(is.na(fg) | fg==0) == 0]

    expect_equal(nrow(snpinfo), ncol(fg))
    expect_equal(rownames(snpinfo), colnames(fg))
})

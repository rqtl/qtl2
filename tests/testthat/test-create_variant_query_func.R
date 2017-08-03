context("create_variant_query_func")

test_that("create_variant_query_func works", {

    # use file name
    dbfile <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2db")
    qf <- create_variant_query_func(dbfile)

    expected <- structure(list(snp_id = c("rs213863525", "rs235573572", "rs253240367"),
                               chr = c("2", "2", "2"), pos = c(97.300098, 97.30014, 97.300197),
                               alleles = c("A|T", "T|C", "G|T"), sdp = c(32L, 64L, 64L),
                               ensembl_gene = c("ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587"),
                               consequence = c("intron_variant", "intron_variant", "intron_variant"),
                               type = c("snp", "snp", "snp")),
                          .Names = c("snp_id", "chr", "pos", "alleles", "sdp", "ensembl_gene", "consequence", "type"),
                          row.names = c(NA, -3L), class = "data.frame")

    expect_equal(qf(2, 97.3, 97.3002), expected)

    # use db connection
    library(RSQLite)
    db <- dbConnect(SQLite(), dbfile)
    qf2 <- create_variant_query_func(db=db)
    expect_equal(qf2(2, 97.3, 97.3002), expected)
    dbDisconnect(db)

    # include filter
    qf3 <- create_variant_query_func(dbfile, filter="snp_id = 'rs235573572'")
    expected_sub <- expected[2,,drop=FALSE]
    rownames(expected_sub) <- 1L
    expect_equal(qf3(2, 0, 200), expected_sub)

})

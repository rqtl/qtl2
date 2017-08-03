context("create_gene_query_func")

test_that("create_gene_query_func works", {

    # use file name
    dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2db")
    qf <- create_gene_query_func(dbfile)

    expected <- structure(list(chr = c("2", "2", "2"),
                               source = c("MGI", "MGI", "MGI"),
                               type = c("gene", "gene", "pseudogene"),
                               start = c(91.407054, 96.318169, 97.947953),
                               stop = c(179.990882, 97.631666, 97.949756),
                               score = c(NA_real_, NA_real_, NA_real_),
                               strand = c("-", "+", "-"),
                               phase = c(NA_real_, NA_real_, NA_real_),
                               ID = c("MGI:MGI:5451848", "MGI:MGI:2442636", "MGI:MGI:3651056"),
                               Name = c("Gm22071", "Lrrc4c", "Gm13803"),
                               Parent = c(NA_character_, NA_character_, NA_character_),
                               Dbxref = c("ENSEMBL:ENSMUSG00000080507", "VEGA:OTTMUSG00000014522,NCBI_Gene:241568,ENSEMBL:ENSMUSG00000050587",
                                          "VEGA:OTTMUSG00000014532,NCBI_Gene:621146,ENSEMBL:ENSMUSG00000082820"),
                               mgiName = c("predicted gene%2c 22071", "leucine rich repeat containing 4C",
                                           "predicted gene 13803"),
                               bioType = c("miRNA gene", "protein coding gene", "pseudogene")),
                          .Names = c("chr", "source", "type", "start", "stop", "score", "strand",
                                     "phase", "ID", "Name", "Parent", "Dbxref", "mgiName", "bioType"),
                          row.names = c(NA, -3L), class = "data.frame")
    expect_equal(qf(2, 97.5, 98.0), expected)

    # use db connection
    library(RSQLite)
    db <- dbConnect(SQLite(), dbfile)
    qf2 <- create_gene_query_func(db=db)
    expect_equal(qf2(2, 97.5, 98.0), expected)
    dbDisconnect(db)

    # include filter
    qf3 <- create_gene_query_func(dbfile, filter="Name = 'Lrrc4c'")
    expected_sub <- expected[2,,drop=FALSE]
    rownames(expected_sub) <- 1L
    expect_equal(qf3(2, 97.5, 98), expected_sub)

})

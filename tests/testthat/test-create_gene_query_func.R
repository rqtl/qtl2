context("create_gene_query_func")

test_that("create_gene_query_func works", {

    # use file name
    dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
    qf <- create_gene_query_func(dbfile)

    expected <- structure(list(chr = c("2", "2", "2"),
                               source = c("MGI", "MGI", "MGI"),
                               type = c("gene", "gene", "pseudogene"),
                               start = c(96.317583, 97.626967, 97.947953),
                               stop = c(97.631672, 97.627052, 97.949756),
                               score = c(NA_real_, NA_real_, NA_real_),
                               strand = c("+", "-", "-"),
                               phase = c(NA_real_, NA_real_, NA_real_),
                               ID = c("MGI_C57BL6J_2442636", "MGI_C57BL6J_5690712", "MGI_C57BL6J_3651056"),
                               Name = c("Lrrc4c", "Gm44320", "Gm13803"),
                               Parent = c(NA_character_, NA_character_, NA_character_),
                               Dbxref = c("ENSEMBL:ENSMUSG00000050587,NCBI_Gene:241568",
                                          "ENSEMBL:ENSMUSG00000105133",
                                          "ENSEMBL:ENSMUSG00000082820,NCBI_Gene:621146"),
                               gene_id = c("MGI:2442636", "MGI:5690712", "MGI:3651056"),
                               mgi_type = c("protein coding gene", "miRNA gene", "pseudogene"),
                               description=c("leucine rich repeat containing 4C",
                                             "predicted gene, 44320",
                                             "predicted gene 13803")),
                          .Names = c("chr", "source", "type", "start", "stop", "score", "strand",
                                     "phase", "ID", "Name", "Parent", "Dbxref", "gene_id", "mgi_type", "description"),
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
    expected_sub <- expected[1,,drop=FALSE]
    rownames(expected_sub) <- 1L
    expect_equal(qf3(2, 97.5, 98), expected_sub)

})

context("top snps from snp association analysis")

test_that("top_snps() works", {
    if(isnt_karl()) skip("this test only run locally")

    # load example DO data from web
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)

    # subset to chr 2
    DOex <- DOex[,"2"]

    # calculate genotype probabilities and convert to allele probabilities
    pr <- calc_genoprob(DOex, error_prob=0.002, cores=0) # multi-core
    apr <- genoprob_to_alleleprob(pr)

    # download snp info from web
    tmpfile <- tempfile()
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/c2_snpinfo.rds")
    download.file(file, tmpfile, quiet=TRUE)
    snpinfo <- readRDS(tmpfile)
    unlink(tmpfile)

    # calculate strain distribution patterns
    snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])

    # identify equivalent SNPs
    snpinfo <- index_snps(DOex$pmap, snpinfo)

    # convert to snp probabilities
    snp_pr <- genoprob_to_snpprob(apr, snpinfo)

    # perform SNP association analysis (here, ignoring residual kinship)
    out_snps <- scan1(snp_pr, DOex$pheno)

    # table with top SNPs
    result <- top_snps(out_snps, snpinfo)

    expected <- structure(list(snp_id = c("rs212414861", "rs229578122", "rs254318131",
                               "rs217679969", "rs238404461", "rs262749925", "rs231282689", "rs260286709",
                               "rs27396282", "rs263841533", "rs231205920", "rs242885221", "rs220351620",
                               "rs52579091", "rs243489710", "rs244316263", "rs219729956", "rs235315566",
                               "rs250167663", "rs234831418", "rs240832432", "rs220815439", "rs579950897",
                               "rs224994851", "rs248208898", "rs245525925", "rs229631954"),
                               chr = c("2", "2", "2", "2", "2", "2", "2", "2", "2", "2",
                               "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2",
                               "2", "2", "2", "2", "2"),
                               pos_Mbp = c(96.754889, 96.754937,
                               96.75494, 96.754955, 96.754963, 96.75497, 96.754979, 96.75518,
                               96.755219, 96.755411, 96.755571, 96.756558, 97.579618, 98.072994,
                               98.110607, 98.125288, 98.125343, 98.204326, 98.214808, 98.217039,
                               98.222529, 98.225114, 98.242395, 98.395104, 98.422488, 98.45411,
                               98.477226),
                               alleles = c("C|G", "T|A", "C|T", "G|T", "T|G",
                               "C|G", "C|G/T", "G|A", "T|C", "T|C", "G|A", "T|C", "A|G",
                               "C|T", "C|T", "C|A", "C|A", "G|T", "C|T", "T|G", "C|A", "G|A",
                               "C|T", "A|T", "C|T", "C|A", "C|T"),
                               AJ =    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               B6 =    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               `129` = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               NOD =   c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               NZO =   c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
                               CAST =  c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               PWK =   c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               WSB =   c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                               sdp = c(48L, 48L, 48L, 48L, 48L, 48L, 48L, 48L, 48L, 48L, 48L, 48L, 16L,
                               16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L, 16L),
                               index = c(3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L, 3264L,
                               16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L, 16182L),
                               interval = rep(c(64L,65L),c(12,15)),
                               on_map=rep(FALSE, 27),
                               lod = c(6.59628910475167, 6.59628910475167, 6.59628910475167, 6.59628910475167, 6.59628910475167,
                               6.59628910475167, 6.59628910475167, 6.59628910475167, 6.59628910475167,
                               6.59628910475167, 6.59628910475167, 6.59628910475167, 5.81902562740654,
                               5.81902562740654, 5.81902562740654, 5.81902562740654, 5.81902562740654,
                               5.81902562740654, 5.81902562740654, 5.81902562740654, 5.81902562740654,
                               5.81902562740654, 5.81902562740654, 5.81902562740654, 5.81902562740654,
                               5.81902562740654, 5.81902562740654)),
                          .Names = c("snp_id",
                          "chr", "pos_Mbp", "alleles", "AJ", "B6", "129", "NOD", "NZO",
                          "CAST", "PWK", "WSB", "sdp", "index", "interval", "on_map", "lod"),
                          row.names = c(3264L,
                          3265L, 3266L, 3267L, 3268L, 3269L, 3271L, 3274L, 3275L, 3283L,
                          3288L, 3289L, 16182L, 22474L, 23017L, 23184L, 23186L, 23360L,
                          23494L, 23518L, 23601L, 23649L, 23830L, 26015L, 26213L, 26554L,
                          26891L), class = "data.frame")

    expect_equal(result, expected)

    # top SNPs among the distinct subset at which calculations were performed
    result <- top_snps(out_snps, snpinfo, show_all_snps=FALSE)
    expect_equal(result, expected[c(1,13),])

    # top SNPs within 1.0 LOD
    result <- top_snps(out_snps, snpinfo, 0.5)
    expect_equal(result, expected[1:12,])

})

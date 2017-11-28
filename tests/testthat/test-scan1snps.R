context("snp association scan")

test_that("scan1snps works", {

    if(isnt_karl()) skip("this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)

    snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
    queryf <- create_variant_query_func(snpdb_file)

    expected <- structure(list(lod = structure(c(0.0974129525672618, 0.910244151499162,
        0.42069693283195, 0.75025701035222, 0.460469686327349, 0.354125542956334,
        0.0073598146805498, 0.0589439809692482, 0.82232343090451, 0.705784666377176,
        0.115892581656003, 0.150541130672339, 0.279942408253708, 0.01042841525841,
        0.187492789850339, 0.331157659586969, 0.0293717397245308, 0.195738520958786,
        0.251698127277011), .Dim = c(19L, 1L), .Dimnames = list(c("rs259104594",
        "rs49962811", "rs250912493", "rs49310503", "rs27363853", "rs27363851",
        "rs33556222", "rs387021772;rs258896118;rs220544684", "rs227317919",
        "rs33673239", "rs33438111", "rs27363804", "rs29764604", "rs241221617",
        "rs27413286", "rs27395539", "rs27395529", "rs27395515", "2:97277447_ATTT/ATT"),
        "OF_immobile_pct"), sample_size = structure(20L, .Names = "OF_immobile_pct"),
        class = c("scan1", "matrix")), snpinfo = structure(list(snp_id = c("rs259104594",
        "rs49962811", "rs250912493", "rs49310503", "rs27363853", "rs27363851",
        "rs33556222", "rs387021772;rs258896118;rs220544684", "rs227317919",
        "rs33673239", "rs33438111", "rs27363804", "rs29764604", "rs241221617",
        "rs27413286", "rs27395539", "rs27395529", "rs27395515", "2:97277447_ATTT/ATT"),
        chr = c("2", "2", "2", "2", "2", "2", "2", "2", "2", "2",
                "2", "2", "2", "2", "2", "2", "2", "2", "2"),
        pos = c(97.20006, 97.200082, 97.200598, 97.200765, 97.2008, 97.200955, 97.201002,
                97.205127, 97.20862, 97.209716, 97.210417, 97.215528, 97.219902,
                97.221925, 97.232949, 97.248752, 97.24921, 97.249953, 97.277447),
        alleles = c("T|A", "T|C", "G|A", "A|C", "T|A", "G|A", "T|G",
                    "CT|CTTTT/CTT/C", "T|C", "C|T", "C|T/G", "C|A/T", "T|C", "A|T",
                    "A|G", "T|C", "T|C", "G|A", "ATT|ATTT/AT/A"),
        sdp = c(32L, 8L, 64L, 96L, 12L, 4L, 108L, 23L, 128L, 17L, 40L, 72L, 2L, 44L, 100L,
                36L, 68L, 76L, 170L),
        ensembl_gene = c("ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000088633", "ENSMUSG00000050587",
                         "ENSMUSG00000050587", "ENSMUSG00000050587", "ENSMUSG00000050587"),
        consequence = c("intron_variant", "intron_variant", "intron_variant",
                        "intron_variant", "intron_variant", "intron_variant", "intron_variant",
                        "intron_variant", "intron_variant", "intron_variant", "intron_variant",
                        "intron_variant", "intron_variant", "intron_variant", "upstream_gene_variant",
                        "intron_variant", "intron_variant", "intron_variant", "intron_variant"),
        A_J = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L, 1L, 2L, 1L, 1L, 1L,
                1L, 1L, 1L, 1L, 1L, 1L),
        C57BL_6J = c(1L, 1L, 1L, 1L, 1L, 1L,
                     1L, 3L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L),
        `129S1_SvImJ` = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L,
                          2L, 1L),
        NOD_ShiLtJ = c(1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L,
                       2L, 2L, 1L, 2L, 1L, 1L, 1L, 2L, 3L),
        NZO_HlLtJ = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 4L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
        CAST_EiJ = c(2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 3L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 4L),
        PWK_PhJ = c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 3L, 1L, 1L, 2L, 1L, 2L, 2L, 1L),
        WSB_EiJ = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L),
        type = c("snp", "snp", "snp", "snp", "snp", "snp", "snp",
                 "indel", "snp", "snp", "snp", "snp", "snp", "snp", "snp", "snp",
                 "snp", "snp", "indel"),
        index = 1:19,
        interval = c(64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L,
                     64L, 64L, 64L, 64L),
        on_map = c(FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE)),
        .Names = c("snp_id", "chr", "pos", "alleles", "sdp", "ensembl_gene", "consequence", "A_J",
                   "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLtJ", "CAST_EiJ",
                   "PWK_PhJ", "WSB_EiJ", "type", "index", "interval", "on_map"),
        row.names = c(1L, 2L, 14L, 22L, 26L, 34L, 36L, 152L, 228L, 256L, 285L, 448L, 543L,
                      569L, 780L, 950L, 963L, 982L, 1324L), class = "data.frame")),
        .Names = c("lod", "snpinfo"))

    # using query function for a defined region
    out <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97.2, end=97.3)
    expect_equal(out, expected)

    # using a pre-defined table of snps
    snpinfo <- queryf(2, 97.2, 97.3)
    out2 <- scan1snps(probs, DOex$pmap, DOex$pheno, snpinfo=snpinfo)
    expect_equal(out, out2)

    # using a pre-defined table of snps
    out3 <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf)
    expect_equal(dim(out3$lod), c(74,1))
    expect_equal(dim(out3$snpinfo), c(74,19))

    # same, keeping all snps
    out3 <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, keep_all_snps=TRUE)
    expect_equal(dim(out3$lod), c(74,1))
    expect_equal(dim(out3$snpinfo), c(13836,19))

})

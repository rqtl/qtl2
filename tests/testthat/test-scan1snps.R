context("snp association scan")

test_that("scan1snps works", {

    skip_if(isnt_karl(), "this test only run locally")

    RNGkind("Mersenne-Twister") # make sure we're using the standard RNG
    set.seed(20180727)

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)

    snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
    queryf <- create_variant_query_func(snpdb_file)

    expected <- list(lod = structure(c(0.0792733888900488, 0.692764338278575, 0.427585764384135,  0.794318319047203,
                                       0.252436101716556,  1.10965528916839,  0.0052242742952302, 0.044886580360961,
                                       0.77913769683386,   0.723589209750504, 0.130713797731086,  0.150391046233995,
                                       0.254799450048782,  0.00724275232762306, 0.19073862233745, 0.27168521752353,
                                       0.0404789019635032, 0.171541072518635, 0.369411627351868),
                                     .Dim = c(19L, 1L),
                                     .Dimnames = list(c("rs259104594", "rs49962811", "rs250912493", "rs49310503",
                                                        "rs27363853", "rs27363851",  "rs33556222",
                                                        "rs387021772;rs258896118;rs220544684",
                                                        "rs227317919", "rs33673239", "rs33438111", "rs27363804",
                                                        "rs29764604", "rs241221617", "rs27413286", "rs27395539",
                                                        "rs27395529", "rs27395515", "2:97277447_ATTT/ATT"),
                                                      "OF_immobile_pct"),
                                     sample_size = c(OF_immobile_pct = 20L), class = c("scan1", "matrix")),
                     snpinfo = structure(list(snp_id = c("rs259104594", "rs49962811", "rs250912493", "rs49310503",
                                                         "rs27363853", "rs27363851", "rs33556222",
                                                         "rs387021772;rs258896118;rs220544684", "rs227317919",
                                                         "rs33673239", "rs33438111", "rs27363804", "rs29764604",
                                                         "rs241221617", "rs27413286", "rs27395539", "rs27395529",
                                                         "rs27395515", "2:97277447_ATTT/ATT"),
                                              chr = rep("2", 19),
                                              pos = c(97.20006, 97.200082, 97.200598, 97.200765, 97.2008,
                                                      97.200955, 97.201002, 97.205127, 97.20862, 97.209716,
                                                      97.210417, 97.215528, 97.219902, 97.221925, 97.232949,
                                                      97.248752, 97.24921, 97.249953, 97.277447),
                                              alleles = c("T|A", "T|C", "G|A", "A|C", "T|A", "G|A", "T|G",
                                                          "CT|CTTTT/CTT/C", "T|C", "C|T", "C|T/G", "C|A/T",
                                                          "T|C", "A|T", "A|G", "T|C", "T|C", "G|A",
                                                          "ATT|ATTT/AT/A"),
                                              sdp = c(32L, 8L, 64L, 96L, 12L, 4L, 108L, 23L, 128L, 17L,
                                                      40L, 72L, 2L, 44L, 100L, 36L, 68L, 76L, 170L),
                                              ensembl_gene = c(rep("ENSMUSG00000050587", 14),
                                                               "ENSMUSG00000088633,ENSMUSG00000050587",
                                                               rep("ENSMUSG00000050587", 4)),
                                              consequence = c(rep("ENSMUSG00000050587:intron_variant", 14),
                                     "ENSMUSG00000088633:upstream_gene_variant,ENSMUSG00000050587:intron_variant",
                                                              rep("ENSMUSG00000050587:intron_variant", 4)),
                                              A_J = c(1, 1, 1, 1, 1, 1, 1, 4, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                              C57BL_6J = c(1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2),
                                              `129S1_SvImJ` = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2,2,2,2,2,1),
                                              NOD_ShiLtJ = c(1, 2, 1, 1, 2,1,2,1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 2, 3),
                                              NZO_HlLtJ = c(1, 1, 1, 1, 1, 1, 1, 4, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                                              CAST_EiJ = c(2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 2, 2, 1, 1, 4),
                                              PWK_PhJ = c(1, 1, 2, 2, 1, 1, 2, 1, 1, 1, 1, 3, 1, 1, 2, 1, 2, 2, 1),
                                              WSB_EiJ = c(1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3),
                                              type = c("snp", "snp", "snp", "snp", "snp", "snp", "snp", "indel",
                                                       "snp", "snp", "snp", "snp", "snp", "snp", "snp", "snp",
                                                       "snp", "snp", "indel"),
                                              index = 1:19,
                                              interval = rep(65L, 19),
                                              on_map = rep(FALSE, 19)),
                                              row.names = c(1L, 2L, 14L, 22L, 26L, 34L, 36L, 152L, 228L, 256L,
                                                            285L, 448L, 543L, 569L, 780L, 950L, 963L, 982L, 1324L),
                                              class = "data.frame"))

    # using query function for a defined region
    out <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97.2, end=97.3)
    expect_equal(out, expected)

    # if probs and map don't conform, should get an error
    junk_map <- DOex$pmap
    junk_map[[1]] <- junk_map[[1]][-1]
    expect_error( scan1snps(probs, junk_map, DOex$pheno, query_func=queryf, chr=2,
                            start=97.2, end=97.3) )

    # using a pre-defined table of snps
    snpinfo <- queryf(2, 97.2, 97.3)
    out2 <- scan1snps(probs, DOex$pmap, DOex$pheno, snpinfo=snpinfo)
    expect_equal(out, out2)

    # using a pre-defined table of snps
    out3 <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf)
    expect_equal(dim(out3$lod), c(49,1))
    expect_equal(dim(out3$snpinfo), c(49,19))

    # same, keeping all snps
    out3 <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, keep_all_snps=TRUE)
    expect_equal(dim(out3$lod), c(49,1))
    expect_equal(dim(out3$snpinfo), c(13836,19))

    # with weights, vs lm()
    w <- setNames(runif(n_ind(DOex), 0, 10), ind_ids(DOex))
    outw <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97.2, end=97.3,
                      weights=w)

    snpinfo <- queryf(2, 97.2, 97.3)
    snpinfo <- index_snps(DOex$pmap, snpinfo)
    snp_probs <- genoprob_to_snpprob(probs, snpinfo)

    p <- pull_genoprobpos(snp_probs, "rs227317919")

    y <- DOex$pheno[rownames(p),,drop=FALSE]
    w <- w[rownames(p)]
    out_lm <- lm(y ~ -1 + p, weights=w)
    out_lm0 <- lm(y ~ 1, weights=w)

    expect_equal(outw$lod["rs227317919", 1], nrow(p)/2*log10(sum(w*out_lm0$resid^2)/sum(w*out_lm$resid^2)))

    # binary trait
    binphe <- setNames((DOex$pheno > quantile(DOex$pheno, 0.8, na.rm=TRUE))*1, rownames(DOex$pheno))
    biny <- binphe[rownames(p)]

    outbin <- scan1snps(probs, DOex$pmap, binphe, query_func=queryf, chr=2, start=97.2, end=97.3,
                        model="binary")

    out_glm <- glm(biny ~ -1 + p, family=binomial(link=logit))
    out_glm0 <- glm(biny ~ 1, family=binomial(link=logit))
    expect_equal(outbin$lod["rs227317919", 1], (out_glm0$deviance - out_glm$deviance)/2/log(10))

    # binary trait with weights
    w <- ceiling(w) # avoid warning from glm() about weights not being integers
    outbinw <- scan1snps(probs, DOex$pmap, binphe, query_func=queryf,
                         chr=2, start=97.2, end=97.3, model="binary",
                         weights=w, eta_max=25, maxit=1000, tol=1e-4)

    out_glm <- glm(biny ~ -1 + p, family=binomial(link=logit), weights=w)
    out_glm0 <- glm(biny ~ 1, family=binomial(link=logit), weights=w)

    expect_equal(outbinw$lod["rs227317919", 1],
                 (out_glm0$deviance - out_glm$deviance)/2/log(10),
                 tolerance=5e-5)

})

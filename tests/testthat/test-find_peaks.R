context("find_peaks lod_int, and bayes_int")

# read data
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

# calculate genotype probabilities
map <- insert_pseudomarkers(iron$gmap, step=1)
probs <- calc_genoprob(iron, map, error_prob=0.002)

# grab phenotypes and covariates; ensure that covariates have names attribute
pheno <- iron$pheno
covar <- match(iron$covar$sex, c("f", "m")) # make numeric
names(covar) <- rownames(iron$covar)
Xcovar <- get_x_covar(iron)

# perform genome scan
out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

test_that("find_peaks works", {

    expected_2_1 <- structure(list(lodindex = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L),
                                   lodcolumn = c("liver", "liver", "liver", "liver", "liver", "liver", "liver",
                                   "liver", "liver", "liver", "liver", "liver", "spleen", "spleen", "spleen"),
                                   chr = structure(c(1L, 2L, 4L, 7L, 7L, 8L, 11L, 13L, 15L, 16L, 17L, 19L, 8L, 9L, 10L),
                                                   .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                                                              "12", "13", "14", "15", "16", "17", "18", "19", "X"), class = "factor"),
                                   pos = c(90.3, 56.8, 10.9, 25.1, 49.1, 41, 26, 32.5, 49.2, 28.6, 4.3, 38.3, 13.6, 56.6, 37),
                                   lod = c(2.53157055941134, 4.85599006729575, 2.74279943572619, 3.77915164598982, 4.42062344751016,
                                   3.49656665028876, 2.61495248604883, 3.40360749696547, 4.03706766736211, 6.35264382891418,
                                   2.71848298281312, 3.0383533142527, 5.35976176735248, 12.5986057120873, 2.23257885452713)),
                              .Names = c("lodindex", "lodcolumn", "chr", "pos", "lod"), row.names = c(NA, 15L),
                              class = "data.frame")

    expect_equal(find_peaks(out, map, 2, 1), expected_2_1)

    # test sort
    expect_equal(find_peaks(out, map, 2, 1, sort_by="pos"), expected_2_1[c(1,2,3,4,5,13,6,14,15,7,8,9,10,11,12),])
    expect_equal(find_peaks(out, map, 2, 1, sort_by="lod"), expected_2_1[c(14,10,13,2,5,9,4,6,8,12,3,11,7,1,15),])

    expected_4_Inf <- expected_2_1[c(2, 5, 9, 10, 13, 14),]
    rownames(expected_4_Inf) <- NULL
    expect_equal(find_peaks(out, map, 4, Inf), expected_4_Inf)

    expected_3_1 <- expected_2_1[c(2,4,5,6,8,9,10,12,13,14),]
    rownames(expected_3_1) <- NULL
    expect_equal(find_peaks(out, map, 3, 1), expected_3_1)

    # separate X threshold
    expect_equal(find_peaks(out, map, 2, 1, thresholdX=1.5),
                 rbind(expected_2_1,
                       data.frame(lodindex=2, lodcolumn="spleen",
                                  chr="X", pos=29.5, lod=1.65321447653698,
                                  stringsAsFactors=FALSE)))

    # column-specific thresholds
    result <- find_peaks(out, map, c(4,2), 1, thresholdX=c(0.3, 1.5), peakdropX=c(0.1, 1))
    expected <- rbind(expected_2_1[c(2,5,9,10),],
                      data.frame(lodindex=rep(1,2), lodcolumn=rep("liver",2),
                                 chr=rep("X",2), pos=c(29.5, 57.9),
                                 lod=c(0.832206643447098, 0.35793505359049),
                                 stringsAsFactors=FALSE),
                      expected_2_1[13:15,],
                      data.frame(lodindex=2, lodcolumn="spleen",
                                 chr="X", pos=29.5, lod=1.65321447653698,
                                 stringsAsFactors=FALSE))
    rownames(expected) <- NULL
    expect_equal(result, expected)

    # lod support intervals
    expect_equal(find_peaks(out, map, 4, 2, 1.5),
                 cbind(expected_4_Inf,
                       ci_lo=c(48.1,  1.1, 16.4,  6.6,  0.0, 43.7),
                       ci_hi=c(73.2, 53.6, 49.2, 40.4, 32.7, 61.2)))

    expect_equal(find_peaks(out, map, 4, 2, 2),
                 cbind(expected_4_Inf,
                       ci_lo=c(48.1,  1.1, 16.4,  6.6,  0.0, 43.7),
                       ci_hi=c(73.2, 53.6, 49.2, 40.4, 32.7, 61.2)))

    expected_3_1_lodint <- cbind(expected_3_1,
                                 ci_lo=c(48.1,  1.1, 37.2, 17.3, 17.5, 16.4,  6.6,  3.3,  0.0, 53.6),
                                 ci_hi=c(73.2, 28.4, 53.6, 69.9, 40.4, 49.2, 40.4, 38.3, 17.3, 61.2))

    expect_equal(find_peaks(out, map, 3, 1, 0.9), expected_3_1_lodint)

    skip_on_cran()

    # sorting
    expect_equal(find_peaks(out, map, 3, 1, 0.9, sort_by="pos"), expected_3_1_lodint[c(1,2,3,9,4,10,5,6,7,8),])
    expect_equal(find_peaks(out, map, 3, 1, 0.9, sort_by="lod"), expected_3_1_lodint[c(10,7,9,1,3,6,2,4,5,8),])

    # expand2markers=FALSE
    expect_equal(find_peaks(out, map, 3, 1, 0.9, expand2markers=FALSE),
                 cbind(expected_3_1,
                       ci_lo=c(53.3, 11.1, 38.1, 24.0, 20.5, 38.4, 21.6, 28.3,  8.0, 53.6),
                       ci_hi=c(66.3, 28.1, 53.6, 53.0, 40.4, 49.2, 32.6, 38.3, 17.0, 59.6)))

    # Bayes intervals
    expect_equal(find_peaks(out, map, 4, 2, prob=0.95),
                 cbind(expected_4_Inf,
                       ci_lo=c(48.1, 13.1, 16.4,  6.6,  0.0, 53.6),
                       ci_hi=c(73.2, 53.6, 49.2, 40.4, 32.7, 61.2)))

    expect_equal(find_peaks(out, map, 4, 2, prob=0.99),
                 cbind(expected_4_Inf,
                       ci_lo=c(48.1,  1.1, 16.4,  6.6,  0.0, 43.7),
                       ci_hi=c(73.2, 53.6, 49.2, 40.4, 32.7, 61.2)))

    expected_3_1_bayesint <- cbind(expected_3_1,
                                   ci_lo=c(48.1, 13.1, 37.2, 32.7, 17.5, 16.4,  6.6,  3.3,  0.0, 53.6),
                                   ci_hi=c(73.2, 28.4, 53.6, 69.9, 40.4, 49.2, 40.4, 38.3, 17.3, 61.2))

    expect_equal(find_peaks(out, map, 3, 1, prob=0.8), expected_3_1_bayesint)

    # sorting
    expect_equal(find_peaks(out, map, 3, 1, prob=0.8, sort_by="pos"), expected_3_1_bayesint[c(1,2,3,9,4,10,5,6,7,8),])
    expect_equal(find_peaks(out, map, 3, 1, prob=0.8, sort_by="lod"), expected_3_1_bayesint[c(10,7,9,1,3,6,2,4,5,8),])

    # expand2markers=FALSE
    expect_equal(find_peaks(out, map, 3, 1, 0.8, expand2markers=FALSE),
                 cbind(expected_3_1,
                       ci_lo=c(53.3, 12.1, 39.1, 26.0, 21.5, 39.4, 25.1, 28.3,  9.0, 53.6),
                       ci_hi=c(65.3, 28.1, 53.6, 52.0, 40.4, 49.2, 32.6, 38.3, 17.0, 59.6)))


    # test that it works if output or map are subsetted
    expect_equal( find_peaks(out, map["2"]), find_peaks(subset(out, map, chr="2"), map) )

    # test that it works if output has rows shuffled
    out_shuffled <- out[sample(1:nrow(out)),,drop=FALSE]
    class(out_shuffled) <- c("scan1", "matrix")
    expect_equal( find_peaks(out_shuffled, map), find_peaks(out, map))

    # test that find_peaks works if there are no peaks above threshold
    blank_output <- structure(list(lodindex = numeric(0),
                                   lodcolumn = character(0),
                                   chr = structure(integer(0), .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"),
                                                   class = "factor"),
                                   pos = numeric(0),
                                   lod = numeric(0)),
                              .Names = c("lodindex", "lodcolumn", "chr", "pos", "lod"),
                              row.names = integer(0), class = "data.frame")

    expect_equal( find_peaks(out, map, threshold=9999), blank_output)

    # test that find_peaks works if there are no peaks above threshold
    # like above, but also requesting LOD or Bayes intervals
    blank_output <- structure(list(lodindex = numeric(0),
                                   lodcolumn = character(0),
                                   chr = structure(integer(0), .Label = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", "X"),
                                                   class = "factor"),
                                   pos = numeric(0),
                                   lod = numeric(0),
                                   ci_lo = numeric(0),
                                   ci_hi = numeric(0)),
                              .Names = c("lodindex", "lodcolumn", "chr", "pos", "lod", "ci_lo", "ci_hi"),
                              row.names = integer(0), class = "data.frame")

    expect_equal( find_peaks(out, map, threshold=9999, drop=2), blank_output)
    expect_equal( find_peaks(out, map, threshold=9999, prob=0.95), blank_output)

})

test_that("lod_int works", {

    expected <- cbind(ci_lo=1.1, pos=49.1, ci_hi=53.6)
    rownames(expected) <- 1
    expect_equal(lod_int(out, map, 7, 1), expected)

    expected[1,c(1,3)] <- c(12.1, 53.6)
    expect_equal(lod_int(out, map, 7, 1, expand2markers=FALSE), expected)

    expected[1,c(1,3)] <- c(37.2, 53.6)
    expect_equal(lod_int(out, map, 7, 1, drop=0.5), expected)

    expected2 <- rbind("1"=c(1.1,25.1,28.4), "2"=c(37.2,49.1,53.6))
    colnames(expected2) <- colnames(expected)
    expect_equal(lod_int(out, map, 7, 1, 0, 1.5, 1), expected2)

    expected3 <- cbind(ci_lo=43.7, pos=56.6, ci_hi=61.2)
    rownames(expected3) <- 1
    expect_equal(lod_int(out, map, 9, 2), expected3)

    expected4 <- cbind(ci_lo=0.0, pos=13.6, ci_hi=32.7)
    rownames(expected4) <- 1
    expect_equal(lod_int(out, map, 8, 2), expected4)

    expect_equal( lod_int(out, map, chr="2"), lod_int(subset(out, map, chr="2"), map, chr="2"))
    expect_equal( lod_int(out, map, chr="2"), lod_int(out, map["2"], chr="2"))

})


test_that("bayes_int works", {

    expected <- cbind(ci_lo=13.1, pos=49.1, ci_hi=53.6)
    rownames(expected) <- 1
    expect_equal(bayes_int(out, map, 7, 1), expected)

    expected[1,c(1,3)] <- c(15.1, 53.6)
    expect_equal(bayes_int(out, map, 7, 1, expand2markers=FALSE), expected)

    expected[1,c(1,3)] <- c(37.2, 53.6)
    expect_equal(bayes_int(out, map, 7, 1, prob=0.8), expected)

    expected2 <- rbind("1"=c(1.1,25.1,28.4), "2"=c(37.2,49.1,53.6))
    colnames(expected2) <- colnames(expected)
    expect_equal(bayes_int(out, map, 7, 1, 0, 1.5, 0.95), expected2)

    expected3 <- cbind(ci_lo=53.6, pos=56.6, ci_hi=61.2)
    rownames(expected3) <- 1
    expect_equal(bayes_int(out, map, 9, 2), expected3)

    expected4 <- cbind(ci_lo=0.0, pos=13.6, ci_hi=32.7)
    rownames(expected4) <- 1
    expect_equal(bayes_int(out, map, 8, 2), expected4)

    expect_equal( bayes_int(out, map, chr="2"), bayes_int(subset(out, map, chr="2"), map, chr="2"))
    expect_equal( bayes_int(out, map, chr="2"), bayes_int(out, map["2"], chr="2"))

})


test_that("lod_int and bayes_int give same results as R/qtl", {

    skip_on_cran()

    out_rqtl <- data.frame(chr=map2chr(map),
                           pos=map2pos(map),
                           unclass(out))
    rownames(out_rqtl) <- map2markernames(map)
    class(out_rqtl) <- c("scanone", "data.frame")

    # 1st lod col: chr 2, 7, 15, 16
    # 2nd lod col: chr 8, 9
    for(lod in 1:2) {
        for(chr in list(c(2,7,15,16), c(8,9))[[lod]]) {

            for(drop in c(1, 1.5, 2)) {
                # expand to markers
                li_rqtl <- qtl::lodint(out_rqtl, lod=lod, chr=chr, expandtomarkers=TRUE, drop=drop)
                li_qtl2 <- lod_int(out, map, lod=lod, chr=chr, drop=drop)
                expect_equivalent(li_rqtl[,2], as.numeric(li_qtl2))

                # don't expand to markers
                li_rqtl <- qtl::lodint(out_rqtl, lod=lod, chr=chr, drop=drop)
                li_qtl2 <- lod_int(out, map, lod=lod, chr=chr, expand2markers=FALSE, drop=drop)
                expect_equivalent(li_rqtl[,2], as.numeric(li_qtl2))
            }

            for(prob in c(0.90, 0.95, 0.99)) {
                # expand to markers
                bi_rqtl <- qtl::bayesint(out_rqtl, lod=lod, chr=chr, expandtomarkers=TRUE, prob=prob)
                bi_qtl2 <- bayes_int(out, map, lod=lod, chr=chr, prob=prob)
               expect_equivalent(bi_rqtl[,2], as.numeric(bi_qtl2))

                # don't expand to markers
                bi_rqtl <- qtl::bayesint(out_rqtl, lod=lod, chr=chr, prob=prob)
                bi_qtl2 <- bayes_int(out, map, lod=lod, chr=chr, expand2markers=FALSE, prob=prob)
                expect_equivalent(bi_rqtl[,2], as.numeric(bi_qtl2))
            }
        }
    }


})

test_that("find_peaks works with snpinfo table", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)

    snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
    queryf <- create_variant_query_func(snpdb_file)

    out <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97.2, end=97.3)

    # test max_scan1()
    expect_equal(max(out$lod, out$snpinfo),
                 data.frame(chr="2",
                            pos=out$snpinfo[6,"pos"],
                            OF_immobile_pct=out$lod[6],
                            row.names=out$snpinfo[6,"snp_id"],
                            stringsAsFactors=FALSE))

    # find_peaks (no peaks above threshold)
    expect_equal(find_peaks(out$lod, out$snpinfo),
                 data.frame(lodindex=numeric(0),
                            lodcolumn=character(0),
                            chr=factor(character(0), "2"),
                            pos=numeric(0),
                            lod=numeric(0),
                            stringsAsFactors=FALSE))
    expect_equal(find_peaks(out$lod, out$snpinfo, drop=1.5),
                 data.frame(lodindex=numeric(0),
                            lodcolumn=character(0),
                            chr=factor(character(0), "2"),
                            pos=numeric(0),
                            lod=numeric(0),
                            ci_lo=numeric(0),
                            ci_hi=numeric(0),
                            stringsAsFactors=FALSE))
    expect_equal(find_peaks(out$lod, out$snpinfo, prob=0.95),
                 data.frame(lodindex=numeric(0),
                            lodcolumn=character(0),
                            chr=factor(character(0), "2"),
                            pos=numeric(0),
                            lod=numeric(0),
                            ci_lo=numeric(0),
                            ci_hi=numeric(0),
                            stringsAsFactors=FALSE))

    # find peaks (one peak above threshold)
    expect_equal(find_peaks(out$lod, out$snpinfo, threshold=0.5),
                 data.frame(lodindex=1,
                            lodcolumn="OF_immobile_pct",
                            chr=factor("2", "2"),
                            pos=out$snpinfo[6,"pos"],
                            lod=out$lod[6],
                            row.names=1L,
                            stringsAsFactors=FALSE))
    expect_equal(find_peaks(out$lod, out$snpinfo, threshold=0.5, drop=1.5),
                 data.frame(lodindex=1,
                            lodcolumn="OF_immobile_pct",
                            chr=factor("2", "2"),
                            pos=out$snpinfo[6,"pos"],
                            lod=out$lod[6],
                            ci_lo=min(out$snpinfo$pos),
                            ci_hi=max(out$snpinfo$pos),
                            row.names=1L,
                            stringsAsFactors=FALSE))
    expect_equal(find_peaks(out$lod, out$snpinfo, threshold=0.5, prob=0.95),
                 data.frame(lodindex=1,
                            lodcolumn="OF_immobile_pct",
                            chr=factor("2", "2"),
                            pos=out$snpinfo[6,"pos"],
                            lod=out$lod[6],
                            ci_lo=min(out$snpinfo$pos),
                            ci_hi=max(out$snpinfo$pos),
                            row.names=1L,
                            stringsAsFactors=FALSE))

})

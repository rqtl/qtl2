context("find_peaks")

test_that("find_peaks works", {

    # load qtl2geno package for data and genoprob calculation
    library(qtl2geno)

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # calculate genotype probabilities
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)

    # grab phenotypes and covariates; ensure that covariates have names attribute
    pheno <- iron$pheno
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # perform genome scan
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    expected_2_1 <- structure(list(lodindex = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L),
                                   lodcolumn = c("liver", "liver", "liver", "liver", "liver", "liver", "liver",
                                   "liver", "liver", "liver", "liver", "liver", "spleen", "spleen", "spleen"),
                                   chr = c("1", "2", "4", "7", "7", "8", "11", "13", "15", "16", "17", "19", "8", "9", "10"),
                                   pos = c(90.3, 56.8, 10.9, 25.1, 49.1, 41, 26, 32.5, 49.2, 28.6, 4.3, 38.3, 13.6, 56.6, 37),
                                   lod = c(2.53157055941134, 4.85599006729575, 2.74279943572619, 3.77915164598982, 4.42062344751016,
                                   3.49656665028876, 2.61495248604883, 3.40360749696547, 4.03706766736211, 6.35264382891418,
                                   2.71848298281312, 3.0383533142527, 5.35976176735248, 12.5986057120873, 2.23257885452713)),
                              .Names = c("lodindex", "lodcolumn", "chr", "pos", "lod"), row.names = c(NA, 15L),
                              class = "data.frame")

    expect_equal(find_peaks(out, 2, 1), expected_2_1)

    expected_4_Inf <- expected_2_1[c(2, 5, 9, 10, 13, 14),]
    rownames(expected_4_Inf) <- NULL
    expect_equal(find_peaks(out, 4, Inf), expected_4_Inf)

    expected_3_1 <- expected_2_1[c(2,4,5,6,8,9,10,12,13,14),]
    rownames(expected_3_1) <- NULL
    expect_equal(find_peaks(out, 3, 1), expected_3_1)

    # separate X threshold
    expect_equal(find_peaks(out, 2, 1, thresholdX=1.5),
                 rbind(expected_2_1,
                       data.frame(lodindex=2, lodcolumn="spleen",
                                  chr="X", pos=29.5, lod=1.65321447653698,
                                  stringsAsFactors=FALSE)))

    # column-specific thresholds
    result <- find_peaks(out, c(4,2), 1, thresholdX=c(0.3, 1.5), peakdropX=c(0.1, 1))
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
    expect_equal(find_peaks(out, 4, 2, 1.5),
                 cbind(expected_4_Inf,
                       ci_lo=c(51.3, 33.1, 34.4, 16.6,  6.0, 51.6),
                       ci_hi=c(69.3, 53.6, 49.2, 33.6, 21.0, 59.6)))

    expect_equal(find_peaks(out, 4, 2, 2),
                 cbind(expected_4_Inf,
                       ci_lo=c(50.3, 29.1, 31.4, 14.6,  4.0, 49.6),
                       ci_hi=c(72.3, 53.6, 49.2, 34.6, 25.0, 60.6)))

    expect_equal(find_peaks(out, 3, 1, 0.9),
                 cbind(expected_3_1,
                       ci_lo=c(53.3, 11.1, 38.1, 24.0, 20.5, 38.4, 21.6, 28.3,  8.0, 53.6),
                       ci_hi=c(66.3, 28.1, 53.6, 53.0, 40.4, 49.2, 32.6, 38.3, 17.0, 59.6)))

})

context("index snps")

test_that("index_snps works", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DO_Recla/recla.zip")
    recla <- read_cross2(file)

    # founder genotypes for a set of SNPs
    snpgeno <- rbind(m1=c(3,1,1,3,1,1,1,1),
                     m2=c(1,3,1,3,1,3,1,3),
                     m3=c(1,1,1,1,3,3,3,3),
                     m4=c(1,3,1,3,1,3,1,3),
                     m5=c(1,1,1,3,1,3,1,1),
                     m6=c(1,1,3,1,1,1,3,1))
    sdp <- calc_sdp(snpgeno)
    snpinfo <- data.frame(chr=c("19", "19", "X", "X", "19", "X"),
                          pos=c(40.36, 40.53, 110.91, 111.21, 56.480843, 56.480843),
                          sdp=sdp,
                          snp=c("m1", "m2", "m3", "m4", "m5", "m6"), stringsAsFactors=FALSE)

    snpinfo_windex <- index_snps(recla$pmap, snpinfo)
    snpinfo19_windex <- index_snps(recla$pmap, snpinfo[snpinfo$chr=="19",])
    snpinfoX_windex <- index_snps(recla$pmap, snpinfo[snpinfo$chr=="X",])

    expected <- cbind(snpinfo[c(1,2,5,6,3,4),],
                      index=c(1:3, 1:3),
                      interval=c(91,91,132,112,237,237),
                      on_map=c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))

    expect_equal(snpinfo_windex, rbind(snpinfo19_windex, snpinfoX_windex))
    expect_equal(snpinfo_windex, expected)
    expect_equal(snpinfo19_windex, expected[1:3,])
    expect_equal(snpinfoX_windex, expected[4:6,])

    # repeat with some redundant snps
    snpinfo_extra <- data.frame(chr=c("19", "19", "X", "X"),
                                pos=c(40.45, 40.37, 56.3, 110.83),
                                sdp=c(170, 9, 68, 170),
                                snp=c("m7", "m8", "m9", "m10"),
                                stringsAsFactors=FALSE)
    rownames(snpinfo_extra) <- snpinfo_extra$snp
    snpinfo <- rbind(snpinfo, snpinfo_extra)

    expected <- cbind(snpinfo[c(1,8,7,2,5,  9,6,10,3,4),],
                      index=c(1,1,3,3,5,  1,1,3,4,3),
                      interval=c(91,91,91,91,132,  112,112,237,237,237),
                      on_map=c(FALSE,FALSE,FALSE,FALSE,TRUE,   FALSE,FALSE,FALSE,FALSE,FALSE))

    snpinfo_windex <- index_snps(recla$pmap, snpinfo)
    snpinfo19_windex <- index_snps(recla$pmap, snpinfo[snpinfo$chr=="19",])
    snpinfoX_windex <- index_snps(recla$pmap, snpinfo[snpinfo$chr=="X",])

    expect_equal(snpinfo_windex, rbind(snpinfo19_windex, snpinfoX_windex))
    expect_equal(snpinfo_windex, expected)
    expect_equal(snpinfo19_windex, expected[1:5,])
    expect_equal(snpinfoX_windex, expected[6:10,])

})

context("reduce to index snps")

test_that("reduce_to_index_snps works", {

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

    snpinfo <- cbind(snpinfo[c(1,2,5,6,3,4),],
                     index=c(1:3, 1:3),
                     interval=c(91,91,132,112,237,237),
                     on_map=c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))

    expect_equal(reduce_to_index_snps(snpinfo), snpinfo)
    expect_equal(reduce_to_index_snps(snpinfo[snpinfo$chr==19,]), snpinfo[snpinfo$chr==19,])
    expect_equal(reduce_to_index_snps(snpinfo[snpinfo$chr=="X",]), snpinfo[snpinfo$chr=="X",])

    # repeat with some redundant snps
    snpinfo_extra <- data.frame(chr=c("19", "19", "X", "X"),
                                pos=c(40.45, 40.37, 56.3, 110.83),
                                sdp=c(170, 9, 68, 170),
                                snp=c("m7", "m8", "m9", "m10"),
                                stringsAsFactors=FALSE)
    rownames(snpinfo_extra) <- snpinfo_extra$snp
    snpinfo_extra <- rbind(snpinfo[,1:4], snpinfo_extra)

    snpinfo_extra <- cbind(snpinfo_extra[c(1,8,7,2,3,  9,4,10,5,6),],
                      index=c(1,1,3,3,5,  1,1,3,4,3),
                      interval=c(91,91,91,91,132,  112,112,237,237,237),
                      on_map=c(FALSE,FALSE,FALSE,FALSE,TRUE,   FALSE,FALSE,FALSE,FALSE,FALSE))

    expected <- snpinfo_extra[c(1,3,5,6,8,9),]
    expected$index <- c(1,2,3,1,2,3)

    expect_equal(reduce_to_index_snps(snpinfo_extra), expected)
    expect_equal(reduce_to_index_snps(snpinfo_extra[snpinfo_extra$chr==19,]), expected[expected$chr==19,])
    expect_equal(reduce_to_index_snps(snpinfo_extra[snpinfo_extra$chr=="X",]), expected[expected$chr=="X",])

})

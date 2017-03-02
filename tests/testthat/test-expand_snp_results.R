context("expand snp association results to include all snps")

test_that("expand_snp_results() works", {

    # this is a contrived example; basically a regression test
    out_snps <- rbind(2.3417662372334, 4.6006054612920, 0.0260449600422, 0.7266312038275, 0.2121385389931)
    dimnames(out_snps) <- list(c("rs221396738", "rs264175039", "rs229722012", "rs227574143", "rs27379133"),
                               "OF_immobile_pct")
    attr(out_snps, "sample_size") <- c(OF_immobile_pct = 261)
    attr(out_snps, "snpinfo") <- list("2"=data.frame(
                                          snp_id=c("rs221396738", "rs264175039", "rs227493750", "rs229722012",
                                                   "rs27379137",  "rs227574143", "rs27379136", "rs216849408",
                                                   "rs240457950", "rs263973399", "rs27379135",  "rs27379134",
                                                   "rs238531859", "rs27379133", "rs235268653", "rs260595137",
                                                   "rs219597400", "rs220959761", "rs241157930", "rs251406434"),
                                          chr=rep("2", 20),
                                          pos_Mbp=c(96.500012, 96.500224, 96.500276, 96.500343, 96.500437,
                                                    96.500668, 96.500749, 96.501004, 96.501091, 96.501096,
                                                    96.501138, 96.501283, 96.501312, 96.501406, 96.501861,
                                                    96.501951, 96.501980, 96.501993, 96.502069, 96.502099),
                                          alleles=c("C|T", "A|C", "C|T", "C|G", "C|T",
                                                    "A|C", "A|C", "A|G", "C|A", "T|A",
                                                    "A|T", "A|C", "A|T", "A|T", "C|T",
                                                    "T|A", "C|T", "G|T", "A|C", "C|G"),
                                          AJ=c(1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1),
                                          B6=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                          "129"=c(1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1),
                                          NOD=c(1,1,1,1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,1,1),
                                          NZO=c(1,3,3,3,1,3,3,3,3,3,3,3,1,3,1,1,1,1,1,1),
                                          CAST=c(3,1,1,3,3,1,3,3,1,1,1,3,3,3,3,3,3,3,3,3),
                                          PWK=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                                          WSB=c(1,1,1,3,1,3,3,3,3,3,3,3,1,3,1,1,1,1,1,1),
                                          sdp=c( 32,  16,  16, 176,  32, 144, 176, 176, 144, 144,
                                                144, 176,  32, 189,  32,  32,  32,  32,  32,  32),
                                          index=c(1, 2, 2, 3, 1, 4, 3, 3, 4, 4, 4, 3, 1, 5, 1, 1, 1, 1, 1, 1),
                                      stringsAsFactors=FALSE) )
    class(out_snps) <- c("scan1", "matrix")

    snp_map <- list("2"=c(rs221396738=96.500012, rs264175039=96.500224, rs229722012=96.500343,
                          rs227574143=96.500668,  rs27379133=96.501406))
    attr(snp_map, "is_x_chr") <- c("2"=FALSE)

    index <- c(1,2,2,3,1,4,3,3,4,4,4,3,1,5,1,1,1,1,1,1)
    expected <- list("lod"=cbind(unclass(out_snps)[,1][index]),
                     "map"=list("2"=attr(out_snps, "snpinfo")[["2"]]$pos_Mbp))
    rownames(expected$lod) <- names(expected$map[["2"]]) <- attr(out_snps, "snpinfo")[["2"]]$snp_id
    colnames(expected$lod) <- "OF_immobile_pct"
    attr(expected$map, "is_x_chr") <- c("2"=FALSE)

    expect_equal(expand_snp_results(out_snps, snp_map), expected)

})

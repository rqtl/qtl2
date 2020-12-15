context("calc raw summaries")

test_that("calculation of raw summaries work", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)
    indIDs <- ind_ids(DOex)
    mnames <- marker_names(DOex)
    founder_ids <- rownames(DOex$founder_geno[[1]])

    g <- do.call("cbind", DOex$geno)
    fg <- do.call("cbind", DOex$founder_geno)
    g[g==0] <- NA
    fg[fg==0] <- NA

    # calc_raw_het()
    het <- calc_raw_het(DOex)
    expect_equal(names(het), indIDs)
    expect_true(all(het >= 0 & het <= 1))
    expect_equal(het, rowMeans(g==2, na.rm=TRUE))

    # calc_raw_maf() by marker
    het <- calc_raw_het(DOex, "marker")
    expect_equal(names(het), mnames)
    expect_true(all(het >= 0 & het <= 1))
    expect_equal(het, colMeans(g==2, na.rm=TRUE))

    # calc_raw_maf()
    maf_ind <- calc_raw_maf(DOex)
    maf_mar <- calc_raw_maf(DOex, "marker")
    expect_equal(names(maf_ind), indIDs)
    expect_equal(names(maf_mar), mnames)
    expect_equal(maf_ind, rowMeans(g==3, na.rm=TRUE) + rowMeans(g==2, na.rm=TRUE)/2)
    expect_equal(maf_mar, colMeans(g==3, na.rm=TRUE) + colMeans(g==2, na.rm=TRUE)/2)

    # calc_raw_founder_maf()
    maf_ind <- calc_raw_founder_maf(DOex)
    maf_mar <- calc_raw_founder_maf(DOex, "marker")
    expect_equal(names(maf_ind), founder_ids)
    expect_equal(names(maf_mar), mnames)
    expect_equal(maf_ind, rowMeans(fg==3, na.rm=TRUE))
    expect_equal(maf_mar, colMeans(fg==3, na.rm=TRUE))

    # calc_raw_geno_freq()
    gf_ind <- calc_raw_geno_freq(DOex)
    gf_mar <- calc_raw_geno_freq(DOex, "marker")
    expect_equal(dimnames(gf_ind), list(indIDs, c("AA", "AB", "BB")))
    expect_equal(dimnames(gf_mar), list(mnames, c("AA", "AB", "BB")))
    expect_equal(gf_ind, cbind(AA=rowMeans(g==1, na.rm=TRUE), AB=rowMeans(g==2, na.rm=TRUE), BB=rowMeans(g==3, na.rm=TRUE)))
    expect_equal(gf_mar, cbind(AA=colMeans(g==1, na.rm=TRUE), AB=colMeans(g==2, na.rm=TRUE), BB=colMeans(g==3, na.rm=TRUE)))

    # multi-core
    gf_ind_2 <- calc_raw_geno_freq(DOex, cores=2)
    gf_mar_2 <- calc_raw_geno_freq(DOex, "marker", cores=2)
    expect_equal(gf_ind_2, gf_ind)
    expect_equal(gf_mar_2, gf_mar)

})

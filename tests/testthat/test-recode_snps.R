context("recode_snps")

test_that("recode_snps works", {

    if(isnt_karl()) skip("This test only run locally")

    # load example data
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/master/DOex/DOex.zip")
    DOex <- read_cross2(file)

    pr <- calc_genoprob(DOex, error_prob=0.002, cores=0)
    nmis_ind <- n_missing(DOex)
    nmis_mar <- n_missing(DOex, "mar")

    DOex <- recode_snps(DOex)

    # no change in number missing
    expect_equal(nmis_ind, n_missing(DOex))
    expect_equal(nmis_mar, n_missing(DOex, "mar"))

    # no change in genotype probabilities
    expect_equal(pr, calc_genoprob(DOex, error_prob=0.002, cores=0))

    fg <- do.call("cbind", DOex$founder_geno)
    fg[fg!=1 & fg!=3] <- NA
    expect_true(all(colMeans(fg==1, na.rm=TRUE) >= 0.5))

})

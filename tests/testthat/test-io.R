context("input/output")

test_that("can read grav2 data", {

    zip_file <- system.file("extdata", "grav2.zip", package="qtl2geno")
    grav2 <- read_cross2(zip_file)

    # check that it contains the same stuff
    expect_equal(sort(names(grav2)),
                 c("alleles", "cross_info", "crosstype", "geno", "gmap", "is_female",
                   "is_x_chr", "pheno", "phenocovar"))

    # check summary
    expected <- structure(list(crosstype = "riself", nlines = 162L, nind = 162L,
                               nchr = 5L, nmar = structure(c(26, 42, 64, 35, 67),
                                          .Names = c("1", "2", "3", "4", "5")),
                               npheno = 241L, ncovar = 0, nphenocovar = 1L,
                               totmar = 234), .Names = c("crosstype", "nlines", "nind",
                                              "nchr", "nmar", "npheno", "ncovar", "nphenocovar", "totmar"),
                          class = c("summary.cross2", "list"))
    expect_equal(summary(grav2), expected)

    # calculate QTL genotype probabilities
    pr <- calc_genoprob(grav2, step=1)

})

test_that("can read iron data", {

    zip_file <- system.file("extdata", "iron.zip", package="qtl2geno")
    iron <- read_cross2(zip_file)

    # check that it contains the same stuff
    expect_equal(sort(names(iron)),
                 c("alleles", "covar", "cross_info", "crosstype", "geno", "gmap",
                   "is_female", "is_x_chr", "pheno", "phenocovar"))

    # check summary
    expected <- structure(list(crosstype = "f2", nlines = 284L, nind = 284L,
                               nchr = 20L, nmar = structure(c(3, 5, 2, 2, 2, 2, 7, 8, 5,
                                           2, 7, 2, 2, 2, 2, 5, 2, 2, 2, 2),
                                           .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                           "11", "12", "13", "14", "15", "16", "17", "18", "19", "X")),
                               npheno = 2L, ncovar = 2L, nphenocovar = 1L, totmar = 66),
                          .Names = c("crosstype", "nlines", "nind", "nchr", "nmar", "npheno", "ncovar",
                          "nphenocovar", "totmar"), class = c("summary.cross2", "list"))
    expect_equal(summary(iron), expected)


    # calculate QTL genotype probabilities
    pr <- calc_genoprob(iron, step=1)

})

context("manipulation of phased genotypes")

test_that("subset.phasedgeno works", {

    skip_if(isnt_karl(), "this test only run locally")

    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    pr <- calc_genoprob(DOex, DOex$gmap, cores=15)
    m <- maxmarg(pr)
    ph <- guess_phase(DOex, m)

    ph_sub <- ph[1:10,]
    expected <- lapply(ph, function(a) a[1:10,,])
    class(expected) <- c("phasedgeno", "list")
    expect_equal(ph_sub, expected)

    ph_sub <- subset.phasedgeno(ph, chr=2)
    expected <- unclass(ph)["2"]
    class(expected) <- c("phasedgeno", "list")
    expect_equal(ph_sub, expected)

    ph_sub <- subset.phasedgeno(ph, ind=11:20, chr=2)
    expected <- lapply(ph, function(a) a[11:20,,])["2"]
    class(expected) <- c("phasedgeno", "list")
    expect_equal(ph_sub, expected)

    # test rbind
    phA <- ph[1:10,]
    phB <- ph[11:20,]
    expect_equal(rbind(phA, phB), ph[1:20,])

    # test cbind
    phA <- ph[,2]
    phB <- ph[,3]
    expect_equal(cbind(phA, phB), ph[,2:3])

})

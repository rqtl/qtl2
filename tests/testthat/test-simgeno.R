context("sim_geno")

# these are really just regression tests

test_that("sim_geno riself", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    RNGkind("Mersenne-Twister")
    set.seed(20150918)

    map <- insert_pseudomarkers(grav2$gmap, step=1)
    dr <- sim_geno(grav2, map, n_draws=2, err=0.002)

    expected <- structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                            1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                            1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                            2L), .Dim = c(135L, 2L),
                          .Dimnames = list(c("PVV4", "c1.loc1",
                          "c1.loc2", "c1.loc3", "c1.loc4", "c1.loc5", "c1.loc6", "AXR-1",
                          "c1.loc7", "c1.loc8", "c1.loc9", "HH.335C-Col/PhyA", "c1.loc10",
                          "c1.loc11", "c1.loc12", "EC.480C", "c1.loc13", "c1.loc14", "c1.loc15",
                          "c1.loc16", "c1.loc17", "c1.loc18", "EC.66C", "c1.loc19", "GD.86L",
                          "c1.loc20", "c1.loc21", "c1.loc22", "c1.loc23", "c1.loc24", "c1.loc25",
                          "CH.160L-Col", "CC.98L-Col/101C", "c1.loc26", "c1.loc27", "c1.loc28",
                          "c1.loc29", "c1.loc30", "c1.loc31", "c1.loc32", "c1.loc33", "c1.loc34",
                          "c1.loc35", "AD.121C", "c1.loc36", "AD.106L-Col", "c1.loc37",
                          "c1.loc38", "c1.loc39", "c1.loc40", "c1.loc41", "c1.loc42", "c1.loc43",
                          "c1.loc44", "GB.112L", "c1.loc45", "c1.loc46", "c1.loc47", "c1.loc48",
                          "c1.loc49", "c1.loc50", "GD.97L", "c1.loc51", "c1.loc52", "c1.loc53",
                          "EG.113L/115C", "c1.loc54", "c1.loc55", "c1.loc56", "c1.loc57",
                          "CD.89C", "c1.loc58", "c1.loc59", "c1.loc60", "c1.loc61", "c1.loc62",
                          "BF.206L-Col", "c1.loc63", "c1.loc64", "c1.loc65", "c1.loc66",
                          "c1.loc67", "CH.200C", "c1.loc68", "c1.loc69", "c1.loc70", "c1.loc71",
                          "c1.loc72", "EC.88C", "c1.loc73", "c1.loc74", "c1.loc75", "c1.loc76",
                          "GD.160C", "c1.loc77", "c1.loc78", "c1.loc79", "c1.loc80", "c1.loc81",
                          "HH.375L", "CH.215L", "c1.loc82", "c1.loc83", "c1.loc84", "c1.loc85",
                          "c1.loc86", "c1.loc87", "c1.loc88", "BF.116C", "c1.loc89", "c1.loc90",
                          "GH.157L-Col", "c1.loc91", "c1.loc92", "c1.loc93", "c1.loc94",
                          "c1.loc95", "c1.loc96", "c1.loc97", "CC.318C", "c1.loc98", "c1.loc99",
                          "CD.173L/175C-Col", "c1.loc100", "c1.loc101", "c1.loc102", "c1.loc103",
                          "c1.loc104", "GH.127L-Col/ADH", "c1.loc105", "c1.loc106", "c1.loc107",
                          "c1.loc108", "c1.loc109", "HH.360L-Col"), NULL))

    expect_equal(dr[[1]][138,,], expected)

})

test_that("sim_geno f2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    RNGkind("Mersenne-Twister")
    set.seed(20150918)
    map <- insert_pseudomarkers(iron$gmap, step=1)
    dr <- sim_geno(iron, map, n_draws=2, err=0.002)

    expected <- structure(c(5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 6L, 4L, 6L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L,
                            4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L, 5L, 4L),
                          .Dim = c(2L, 30L, 2L),
                          .Dimnames = list(c("145", "146"), c("DXMit16", "cX.loc30.5",
                          "cX.loc31.5", "cX.loc32.5", "cX.loc33.5", "cX.loc34.5", "cX.loc35.5", "cX.loc36.5",
                          "cX.loc37.5", "cX.loc38.5", "cX.loc39.5", "cX.loc40.5", "cX.loc41.5", "cX.loc42.5",
                          "cX.loc43.5", "cX.loc44.5", "cX.loc45.5", "cX.loc46.5", "cX.loc47.5", "cX.loc48.5",
                          "cX.loc49.5", "cX.loc50.5", "cX.loc51.5", "cX.loc52.5", "cX.loc53.5", "cX.loc54.5",
                          "cX.loc55.5", "cX.loc56.5", "cX.loc57.5", "DXMit186"), NULL))

    expect_equal(dr[["X"]][145:146,,], expected)

})

test_that("sim_geno works when multi-core", {
    if(isnt_karl()) skip("this test only run locally")

    # can't really tell if I'm getting the same answers as w/o multi-core
    # but I can tell if I get the same thing twise, when run from same seed
    # (need to use RNGkind)

    library(qtl)
    data(hyper)
    hyper2 <- convert2cross2(hyper)

    RNGkind("L'Ecuyer-CMRG")
    set.seed(20151209)
    dr <- sim_geno(hyper2, n_draws=2, error_prob=0.002, cores=4)
    set.seed(20151209)
    dr_mc <- sim_geno(hyper2, n_draws=2, error_prob=0.002, cores=4)
    expect_equal(dr_mc, dr)

    data(listeria)
    listeria2 <- convert2cross2(listeria)
    set.seed(20151209)
    dr <- sim_geno(listeria2, n_draws=2, error_prob=0.002, cores=4)
    set.seed(20151209)
    dr_mc <- sim_geno(listeria2, n_draws=2, error_prob=0.002, cores=4)
    expect_equal(dr_mc, dr)

    # re-set RNGkind
    RNGkind("Mersenne-Twister")
})


test_that("sim_geno riself gives same result for lowmem=TRUE and =FALSE", {

    RNGkind("Mersenne-Twister")
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    set.seed(20150918)
    dr <- sim_geno(grav2, map, n_draws=2, err=0.002, lowmem=FALSE)
    set.seed(20150918)
    dr2 <- sim_geno(grav2, map, n_draws=2, err=0.002, lowmem=TRUE)

    expect_equal(dr, dr2)

})

test_that("sim_geno f2 gives same result for lowmem=TRUE and =FALSE", {

    RNGkind("Mersenne-Twister")
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    # order individuals to be sure that the get run in the same order
    p <- paste(iron$is_female, apply(iron$cross_info, 1, paste, collapse=":"), sep=":")
    iron <- iron[order(p),]

    map <- insert_pseudomarkers(iron$gmap, step=1)

    set.seed(20150918)
    dr <- sim_geno(iron, map, n_draws=2, err=0.002, lowmem=FALSE)
    set.seed(20150918)
    dr2 <- sim_geno(iron, map, n_draws=2, err=0.002, lowmem=TRUE)

    expect_equal(dr, dr2)

})

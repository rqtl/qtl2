context("n_missing and n_typed")

test_that("n_missing and n_typed work for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    expect_equal(n_missing(iron[1:4,]), c("1"=0,  "2"=36, "3"=0,  "4"=0))
    expect_equal(n_typed(iron[1:4,]),   c("1"=66, "2"=30, "3"=66, "4"=66))
    expect_equal(n_missing(iron[1:4,], sum="prop"), c("1"=0, "2"=18/33, "3"=0, "4"=0))
    expect_equal(n_typed(iron[1:4,], sum="prop"),   c("1"=1, "2"=15/33, "3"=1, "4"=1))

    expect_equal(n_missing(iron[,3], "mar"),         c("D3Mit22"=129,"D3Mit18"=129))
    expect_equal(n_typed(iron[,3], "mar"),           c("D3Mit22"=155,"D3Mit18"=155))
    expect_equal(n_missing(iron[,3], "mar", "prop"), c("D3Mit22"=129/284,"D3Mit18"=129/284))
    expect_equal(n_typed(iron[,3], "mar", "prop"),   c("D3Mit22"=155/284,"D3Mit18"=155/284))

})

test_that("n_missing and n_typed work for grav2", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))

    expect_equal(n_missing(grav2[125:128,]), c("125"=2,  "126"=3, "127"=4,  "128"=1))
    expect_equal(n_typed(grav2[125:128,]),   c("125"=232,  "126"=231, "127"=230,  "128"=233))
    expect_equal(n_missing(grav2[125:128,], sum="prop"), c("125"=2/234, "126"=3/234, "127"=4/234, "128"=1/234))
    expect_equal(n_typed(grav2[125:128,], sum="prop"), c("125"=232/234, "126"=231/234, "127"=230/234, "128"=233/234))

    expected <- structure(c(1, 2, 2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                            3, 0, 0, 1, 1, 0, 1, 0, 4, 0, 0, 0, 1, 0, 3, 0, 2, 3, 0, 1, 2,
                            1, 4, 2, 3), .Names = c("AD.156C", "BF.325L", "GH.580L", "DF.225L",
                                         "AD.77L", "CH.266C", "CH.610C", "HH.258C", "BH.145C", "BF.226C/BH.58L",
                                         "FD.226C", "GD.145C", "GH.94L", "BF.82C", "GD.465C", "FD.306L",
                                         "EC.495C-Col", "BH.460L", "FD.81L", "BF.105C", "CH.284C", "FD.222L-Col",
                                         "CD.245L", "EG.66L", "CH.65C", "CH.1500C", "BF.221L", "FD.85C",
                                         "GB.150L-Col", "FD.150C", "GD.460L-Col", "CC.332C", "Erecta",
                                         "CH.145L-Col/150C", "AD.191L-Col", "BH.195L-Col", "GD.298C",
                                         "GH.247L", "BH.120L-Col", "DF.140C", "EG.357C/359L-Col",
                                         "EC.235L-Col/247C"))
    expect_equal(n_missing(grav2[,2], "mar"), expected)

})

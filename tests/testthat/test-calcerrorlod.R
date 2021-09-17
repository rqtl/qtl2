context("calc_errorlod")

test_that("calc_errorlod works for an intercross", {

    # this is not much more than a regression test, really
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    pr <- calc_genoprob(iron, err=0.002, map_function="c-f")
    err <- calc_errorlod(iron, pr)

    expected_1to5_chr1 <- structure(c(-2.82278691244427, 0, -2.82278691244427, -3.24248113822165,
                                      0, -2.76132216098372, 0, -2.76132216098372, -3.23746097636371,
                                      0, -3.00415918819379, 0, -3.00415918819379, -2.99576369311289,
                                      0), .Dim = c(5L, 3L),
                                    .Dimnames = list(c("1", "2", "3", "4", "5"),
                                                     c("D1Mit18", "D1Mit80", "D1Mit17")))
    expect_equal(err[[1]][1:5,], expected_1to5_chr1)

    ind <- c(1:5, 244:249)
    expected_chrX <- structure(c(-3.10992625978842, 0, -3.10992625978842, -3.10992625978842,
                                 0, -3.10992625978842, -3.10992625978842, 0, 0, -3.10992625978842,
                                 -3.10992625978842, -3.10992625978842, 0, -3.10992625978842, -3.10992625978842,
                                 0, -3.10992625978842, -3.10992625978842, 0, 0, -3.10992625978842,
                                 -3.10992625978842), .Dim = c(11L, 2L),
                               .Dimnames = list(c("1", "2", "3", "4", "5", "244", "245", "246", "247", "248", "249"),
                                                c("DXMit16", "DXMit186")))
    expect_equal(err$X[ind,], expected_chrX)

    # missing values -> error lod == 0
    for(i in seq(along=err)) {
        if(any(iron$geno[[i]]==0)) {
            expect_true(max(abs(err[[i]][iron$geno[[i]]==0])) < 1e-15)
        }
    }

    # check a locus on chr 3
    p <- pr[[7]][,,4]
    g <- iron$geno[[7]][,4]
    expected <- g
    expected[g==1] <- log10((1-p[,1])/p[,1]*0.25/0.75)[g==1]
    expected[g==2] <- log10((1-p[,2])/p[,2])[g==2]
    expected[g==3] <- log10((1-p[,3])/p[,3]*0.25/0.75)[g==3]

    expect_equal(err[[7]][,4], expected)

    # check a locus on the X chr
    p <- pr$X[,,2]
    g <- iron$geno$X[,2]
    male <- (iron$covar$sex=="m") # genotypes 1/3
    female_f <- (iron$covar$sex=="f" & iron$covar$cross_direction=="(SxB)x(SxB)") # genotypes 1/2
    female_r <- (iron$covar$sex=="f" & iron$covar$cross_direction=="(BxS)x(BxS)") # genotypes 2/3
    expected <- g
    expected[g==1 & male] <- log10((1-p[,5])/p[,5])[g==1 & male]
    expected[g==3 & male] <- log10((1-p[,6])/p[,6])[g==3 & male]
    expected[g==1 & female_f] <- log10((1-p[,1])/p[,1])[g==1 & female_f]
    expected[g==2 & female_f] <- log10((1-p[,2])/p[,2])[g==2 & female_f]
    expected[g==2 & female_r] <- log10((1-p[,3])/p[,3])[g==2 & female_r]
    expected[g==3 & female_r] <- log10((1-p[,4])/p[,4])[g==3 & female_r]

    expect_equal(err$X[,2], expected)

})


test_that("calc_errorlod works for RIL", {

    # this is not much more than a regression test, really
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    pr <- calc_genoprob(grav2, err=0.002, map_function="c-f")
    err <- calc_errorlod(grav2, pr)

    ind <- c(9L, 23L, 25L, 42L, 44L, 49L, 54L, 66L, 134L, 158L)
    mar <- c(32L, 35L, 36L, 48L, 56L)

    expected <- structure(c(-5.82983567140969, -5.82983567140969, -5.82983567140969,
                            -5.82983567140969, -5.82983567140969, -5.82983288817781, -5.82983567140969,
                            -5.82983567140969, -5.82983567140043, -5.82983567140969, -5.60155654532719,
                            -5.60155654532719, -5.60155654532719, -5.60155654532719, -5.60155654532719,
                            -2.61600519342757, -5.60155654532719, -5.60155654532719, -5.60016580984213,
                            -5.60155654266195, -8.29839428989781, -8.29839428989781, -8.29839428989781,
                            -8.29839428989781, -8.29839428989781, -5.31493687487081, -8.29839428989475,
                            -8.29839428989775, -5.47724186085683, -8.29784500171819, -14.3504722627229,
                            -14.350472262709, -14.3504722627229, -14.3504722627229, -14.3504722627229,
                            -14.3504722627229, -14.3504722578884, -14.3504499660042, -14.3504722627229,
                            -14.3504722627119, -5.86578939924429, -4.92654601370327, -5.86578939924429,
                            -7.62403763042528, -7.62403763042528, -5.86578939924429, -7.62403763042528,
                            -7.62392376500139, -7.62403763042527, -7.62392376500139),
                          .Dim = c(10L, 5L),
                          .Dimnames = list(c("9", "23", "25", "42", "44", "49", "54", "66", "134", "158"),
                                           c("DF.300C", "CD.179L", "CH.88L", "BH.81L-Col", "FD.345C")))

    expect_equal(err[[5]][ind,mar], expected)

    # missing values -> error lod == 0
    for(i in seq(along=err)) {
        if(any(grav2$geno[[i]]==0)) {
            expect_true(max(abs(err[[i]][grav2$geno[[i]]==0])) < 1e-15)
        }
    }

    # check a locus on chr 4
    p <- pr[[4]][,,35]
    g <- grav2$geno[[4]][,35]
    expected <- g
    for(i in 1:2) expected[g==i] <- log10((1-p[,i])/p[,i])[g==i]
    expect_equal(err[[4]][,35], expected)

})


test_that("calc_errorlod works for RIL", {

    skip_if(isnt_karl(), "this test only run locally")

    # this is not much more than a regression test, really
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DO_Recla/recla.zip")
    recla <- read_cross2(file)
    recla <- recla[c(1:2,53:54), c("19","X")] # subset to 4 mice and 2 chromosomes
    probs <- calc_genoprob(recla, err=0.002)

    err <- calc_errorlod(recla, probs)

    expected_c19 <- structure(c(-9.32903376530038, -8.47970241856038, 2.32305910951409,
                                -5.89845329201031, -9.47388004360523, -8.33322922086628, -3.99257600876601,
                                -4.20445567758459, -7.50493914753867, -9.14775695748182, -9.09524930864587,
                                -6.79076628467775, -5.15284031339008, -6.99865791362416, -9.76892294431087,
                                -5.40229996503646), .Dim = c(4L, 4L),
                              .Dimnames = list(c("1", "4", "58", "59"),
                                               c("UNC190139327", "UNC190279660", "backupUNC190139992", "UNC190140427")))

    expect_equal(err[[1]][,7:10], expected_c19)

    expected_cX <- structure(c(-6.08215766199779, -18.6736643714498, -6.5183772856749,
                               -7.07293927915729, 0, -16.5645778105845, -0.442200926998983,
                               -1.14542766590646, -17.9466342307662, -12.3112920811831, -16.8345767687497,
                               7.53542542841645, -16.9579570002804, -16.9636066624958, -18.3869460431974,
                               -9.65729043655822, -18.2175770053576, -9.87948829125263, -17.2496166243865,
                               -7.59183451077769), .Dim = 4:5, .Dimnames = list(c("1", "4",
                               "58", "59"), c("UNC200034449", "UNC200261375", "UNC200062736",
                               "backupUNC200189557", "UNC200192385")))

    expect_equal(err$X[,266:270], expected_cX)

})

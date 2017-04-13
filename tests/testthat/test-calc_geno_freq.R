context("Calculate genotype frequencies")

test_that("calc_geno_freq works for an intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    pr <- calc_genoprob(iron, err=0.002)

    # omit_x works?
    expect_equal(calc_geno_freq(pr[,1:19]), calc_geno_freq(pr, omit_x=TRUE))
    expect_equal(calc_geno_freq(pr[,1:19], "marker"), calc_geno_freq(pr, "marker", omit_x=TRUE))

    # by individual
    expected <- structure(list(A = structure(c(0.265431865621607, 0.179999497125949,
                                               0.437768300719857, 0.530653642640462,
                                               0.296799833658536, 0.289346860233589),
                                             .Dim = 2:3, .Dimnames = list(c("88", "105"), c("SS", "SB", "BB"))),
                               X = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5),
                                             .Dim = c(2L, 6L), .Dimnames = list(c("88", "105"),
                                                                                c("SS", "SB", "BS", "BB",
                                                                                  "SY", "BY")))),
                          .Names = c("A", "X"))
    expect_equal(calc_geno_freq(pr[c(88,105),]), expected$A)
    expect_equal(calc_geno_freq(pr[c(88,105),], omit_x=TRUE), expected$A)
    expect_equal(calc_geno_freq(pr[c(88,105),], omit_x=FALSE), expected)

    # by marker
    expected <- structure(list(A = structure(c(0.268363278479528, 0.278689258334566,
                                               0.48461632550778, 0.477608930381241,
                                               0.247020396012691, 0.243701811284194),
                                             .Dim = 2:3, .Dimnames = list(c("D19Mit68", "D19Mit37"),
                                                                          c("SS", "SB", "BB"))),
                               X = structure(c(0.121465495590373, 0.117994689355019,
                                               0.139097884691317, 0.142568690926671,
                                               0.112701216600633, 0.116172022835987,
                                               0.116172022835987, 0.112701216600633,
                                               0.239547079613345, 0.253430304554763,
                                               0.271016300668345, 0.257133075726927),
                                             .Dim = c(2L, 6L), .Dimnames = list(c("DXMit16", "DXMit186"),
                                                                                c("SS", "SB", "BS", "BB",
                                                                                  "SY", "BY")))),
                          .Names = c("A", "X"))
    expect_equal(calc_geno_freq(pr[,c(19,"X")], "marker"), expected$A)
    expect_equal(calc_geno_freq(pr[,c(19,"X")], "marker", omit_x=TRUE), expected$A)
    expect_equal(calc_geno_freq(pr[,c(19,"X")], "marker", omit_x=FALSE), expected)

})

context("Calculate heterozygosities")

test_that("calc_het works for an intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    pr <- calc_genoprob(iron, err=0.002)

    # omit_x works?
    expect_equal(calc_het(pr[,1:19]), calc_het(pr, omit_x=TRUE))
    expect_equal(calc_het(pr[,1:19], "marker"), calc_het(pr, "marker", omit_x=TRUE))

    # by individual
    ind <- c(88,105,231)
    expected <- c(`88` = 0.437768300719857, `105` = 0.530653642640462, `231` = 0.51155729816745)
    expect_equal(calc_het(pr[ind,]), expected)
    expect_equal(calc_het(pr[ind,], omit_x=TRUE), expected)
    expected_wX <- c(`88` = 0.424502594637437, `105` = 0.514573229227115, `231` = 0.526341691848004)
    expect_equal(calc_het(pr[ind,], omit_x=FALSE), expected_wX)

    # by marker
    expected <- c(D19Mit68 = 0.48461632550778, D19Mit37 = 0.477608930381241,
                  DXMit16 = 0.251799101291949, DXMit186 = 0.258740713762659)
    expect_equal(calc_het(pr[,c(19,"X")], "marker"), expected[1:2])
    expect_equal(calc_het(pr[,c(19,"X")], "marker", omit_x=TRUE), expected[1:2])
    expect_equal(calc_het(pr[,c(19,"X")], "marker", omit_x=FALSE), expected)

})

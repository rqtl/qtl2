context("count_xo and locate_xo")

test_that("count_xo and locate_xo work for intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[1:5, c(7,8,"X")]
    pr <- calc_genoprob(iron, error_prob=0.002, map_function="c-f")
    v <- maxmarg(pr)

    n_xo <- count_xo(v)
    expected <- structure(c(2L, 0L, 3L, 2L, 0L, 1L, 0L, 2L, 0L, 1L, 0L, 0L, 0L,
                            0L, 0L), .Dim = c(5L, 3L),
                          .Dimnames = list(c("1", "2", "3", "4", "5"), c("7", "8", "X")))
    expect_equal(n_xo, expected)

    pos_xo <- locate_xo(v, iron$gmap)
    pos_xoB <- locate_xo(v, iron$gmap[c("7", "8", "X")])
    expect_equal(pos_xo, pos_xoB)

    midpts <- sapply(iron$gmap[c(7,8,"X")], function(a) setNames((a[-1] + a[-length(a)])/2, NULL))
    expected <- list("7"=list("1"=midpts[[1]][c(4,6)],
                              "2"=numeric(0),
                              "3"=midpts[[1]][c(3,3,6)],
                              "4"=midpts[[1]][c(3,6)],
                              "5"=numeric(0)),
                     "8"=list("1"=midpts[[2]][3],
                              "2"=numeric(0),
                              "3"=midpts[[2]][c(4,4)],
                              "4"=numeric(0),
                              "5"=midpts[[2]][6]),
                     "X"=list("1"=numeric(0),
                              "2"=numeric(0),
                              "3"=numeric(0),
                              "4"=numeric(0),
                              "5"=numeric(0)))
    expect_equal(pos_xo, expected)

    n_xo_derived <- sapply(pos_xo, sapply, length)
    expect_equal(n_xo, n_xo_derived)


    # sim_geno
    set.seed(85309395)
    dr <- sim_geno(iron, n_draws=4, error_prob=0.002, map_function="c-f")
    n_xo <- count_xo(dr)
    expect_equal(dim(n_xo), c(5,3,4))

    expected <- structure(c(2L, 2L, 3L, 2L, 4L, 3L, 0L, 2L, 0L, 3L, 0L, 0L, 0L,
                            0L, 0L, 2L, 3L, 3L, 2L, 1L, 3L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L,
                            0L, 2L, 0L, 3L, 2L, 0L, 3L, 0L, 2L, 0L, 3L, 0L, 1L, 0L, 0L, 1L,
                            2L, 2L, 3L, 2L, 1L, 2L, 0L, 2L, 0L, 3L, 0L, 0L, 0L, 0L, 1L),
                          .Dim = c(5L, 3L, 4L), .Dimnames = list(c("1", "2", "3", "4", "5"),
                                                                 c("7", "8", "X"), NULL))
    expect_equal(n_xo, expected)

})

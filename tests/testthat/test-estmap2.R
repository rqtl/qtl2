context("est_map gives same results w/ lowmem=T or F")

test_that("est_map2 works with intercross", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(2,3,19,"X")]

    map_lomem <- est_map(iron, lowmem=TRUE)
    map_himem <- est_map(iron, lowmem=FALSE)
    expect_equal(map_himem, map_lomem)

})

test_that("est_map2 works with RIL by selfing", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[,4:5]

    map_lomem <- est_map(grav2, lowmem=TRUE)
    map_himem <- est_map(grav2, lowmem=FALSE)
    expect_equal(map_himem, map_lomem)

})


test_that("est_map2 works with a backcross with both sexes", {

    library(qtl)
    data(hyper)
    set.seed(35288832)
    hyper$pheno$sex <- factor(sample(c("female", "male"), nind(hyper), replace=TRUE))

    hyper <- convert2cross2(hyper)
    hyper <- hyper[,c(4,19,"X")]


    map_lomem <- est_map(hyper, lowmem=TRUE)
    map_himem <- est_map(hyper, lowmem=FALSE)
    expect_equal(map_himem, map_lomem)

})

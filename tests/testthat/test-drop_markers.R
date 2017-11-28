context("drop_markers, drop_nullmarkers, pull_markers")

test_that("drop_markers works with riself", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    # extra markers in list
    markers2drop <- c("BH.342C/347L-Col", "GH.94L", "EG.357C/359L-Col", "blah", "CD.245L", "karl", "caleb", "ANL2")
    expect_warning(drop_markers(grav2, markers2drop), "Some markers not found: blah, karl, caleb")

    # drop some markers
    markers2drop <- markers2drop[-c(4,6,7)]
    result <- drop_markers(grav2, markers2drop)
    expected <- grav2
    for(i in c(2,4)) {
        mn <- marker_names(grav2[,i])
        expected$geno[[i]] <- expected$geno[[i]][,!(mn %in% markers2drop)]
        expected$gmap[[i]] <- expected$gmap[[i]][!(mn %in% markers2drop)]
    }
    expect_equal(result, expected)

    # drop a whole chromosome
    mn <- sample(marker_names(grav2[,3]))
    expect_equal(drop_markers(grav2, mn), grav2[,-3])

})

test_that("drop_markers works with f2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    # drop some markers
    markers2drop <- c("D1Mit80", "D8Mit195", "D3Mit18", "D8Mit31", "D8Mit124",
                      "D5Mit30", "D9Mit10", "D4Mit352", "D1Mit17", "D16Mit70")
    result <- drop_markers(iron, markers2drop)
    expected <- iron
    for(i in c(1,3,4,5,8,9,16)) {
        mn <- marker_names(iron[,i])
        expected$geno[[i]] <- expected$geno[[i]][,!(mn %in% markers2drop),drop=FALSE]
        expected$gmap[[i]] <- expected$gmap[[i]][!(mn %in% markers2drop)]
        expected$pmap[[i]] <- expected$pmap[[i]][!(mn %in% markers2drop)]
    }
    expect_equal(result, expected)

    # drop a few chromosomes
    chr2drop <- c(5,6,9,13)
    mn <- sample(marker_names(iron[,chr2drop]))
    expect_equal(drop_markers(iron, mn), iron[,-chr2drop])

})

test_that("drop_nullmarkers works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    iron$geno[[1]][,1:3] <- 0
    iron$geno[[17]][,1] <- 0

    mn <- c(marker_names(iron[,1]),
            marker_names(iron[,17])[1])

    expect_equal(drop_nullmarkers(iron), drop_markers(iron, mn))

})

test_that("pull_markers works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

    mn <- marker_names(iron)
    to_keep <- sample(mn, 50)
    to_drop <- mn[!(mn %in% to_keep)]

    expect_equal(pull_markers(iron, to_keep), drop_markers(iron, to_drop))

})

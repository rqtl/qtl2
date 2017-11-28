context("rbind and cbind viterbi")

test_that("rbind.viterbi works for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    g <- viterbi(grav2, map, error_prob=0.002)
    gA <- g[1:5,]
    gB <- g[6:12,]
    gC <- g[13:20,]
    gAB <- g[1:12,]
    gABC <- g[1:20,]
    gBACA <- g[c(6:12, 1:5, 13:20, 1:5),]

    expect_equal(rbind(gA, gB), gAB)
    expect_equal(rbind(gA, gB, gC), gABC)
    expect_equal(rbind(gB, gA, gC, gA), gBACA)

})


test_that("rbind.viterbi works for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    g <- viterbi(iron, map, error_prob=0.002)
    gA <- g[2:20,]
    gB <- g[41:60,]
    gC <- g[102:201,]
    gAB <- g[c(2:20,41:60),]
    gABC <- g[c(2:20,41:60,102:201),]
    gBACA <- g[c(41:60,2:20,102:201,2:20),]

    expect_equal(rbind(gA, gB), gAB)
    expect_equal(rbind(gA, gB, gC), gABC)
    expect_equal(rbind(gB, gA, gC, gA), gBACA)

})

test_that("cbind.viterbi works for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    g <- viterbi(grav2[1:10,], map, error_prob=0.002)
    gA <- g[,1:2]
    gB <- g[,5]
    gC <- g[,3:4]
    gAB <- g[,c(1:2,5)]
    gABC <- g[,c(1,2,5,3,4)]
    gBACA <- g[,c(5,1:2,3:4,1:2)]

    expect_equal(cbind(gA, gB), gAB)
    expect_equal(cbind(gA, gB, gC), gABC)
    expect_equal(cbind(gB, gA, gC, gA), gBACA)

})


test_that("cbind.viterbi works for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    g <- viterbi(iron[6:21,], map, error_prob=0.002)
    gA <- g[,2:3]
    gB <- g[,c(4,5,8)]
    gC <- g[,c(19,"X")]
    gAB <- g[,c(2:3,4,5,8)]
    gABC <- g[,c(2:3,4,5,8,19,"X")]
    gBACA <- g[,c(4,5,8,2,3,19,"X",2,3)]

    expect_equal(cbind(gA, gB), gAB)
    expect_equal(cbind(gA, gB, gC), gABC)
    expect_equal(cbind(gB, gA, gC, gA), gBACA)

})

context("rbind and cbind sim_geno")

test_that("rbind.sim_geno works for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[1:30,1:2]
    map <- insert_pseudomarkers(grav2$gmap, step=5)
    draws <- sim_geno(grav2, map, error_prob=0.002, n_draws=8)
    drawsA <- draws[1:5,]
    drawsB <- draws[6:12,]
    drawsC <- draws[13:20,]
    drawsAB <- draws[1:12,]
    drawsABC <- draws[1:20,]
    drawsBACA <- draws[c(6:12, 1:5, 13:20, 1:5),]

    expect_equal(rbind(drawsA, drawsB), drawsAB)
    expect_equal(rbind(drawsA, drawsB, drawsC), drawsABC)
    expect_equal(rbind(drawsB, drawsA, drawsC, drawsA), drawsBACA)

})


test_that("rbind.sim_geno works for iron", {

    skip_on_cran()

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    draws <- sim_geno(iron, map, error_prob=0.002, n_draws=5)
    drawsA <- draws[2:20,]
    drawsB <- draws[41:60,]
    drawsC <- draws[102:201,]
    drawsAB <- draws[c(2:20,41:60),]
    drawsABC <- draws[c(2:20,41:60,102:201),]
    drawsBACA <- draws[c(41:60,2:20,102:201,2:20),]

    expect_equal(rbind(drawsA, drawsB), drawsAB)
    expect_equal(rbind(drawsA, drawsB, drawsC), drawsABC)
    expect_equal(rbind(drawsB, drawsA, drawsC, drawsA), drawsBACA)

})

test_that("cbind.simgeno for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[1:30,]
    map <- insert_pseudomarkers(grav2$gmap, step=5)
    draws <- sim_geno(grav2[1:10,], map, error_prob=0.002, n_draws=8)
    drawsA <- draws[,1:2]
    drawsB <- draws[,5]
    drawsC <- draws[,3:4]
    drawsAB <- draws[,c(1:2,5)]
    drawsABC <- draws[,c(1,2,5,3,4)]
    drawsBACA <- draws[,c(5,1:2,3:4,1:2)]

    expect_equal(cbind(drawsA, drawsB), drawsAB)
    expect_equal(cbind(drawsA, drawsB, drawsC), drawsABC)
    expect_equal(cbind(drawsB, drawsA, drawsC, drawsA), drawsBACA)

})


test_that("cbind.sim_geno works for iron", {

    skip_on_cran()

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    draws <- sim_geno(iron[6:21,], map, error_prob=0.002, n_draws=6)
    drawsA <- draws[,2:3]
    drawsB <- draws[,c(4,5,8)]
    drawsC <- draws[,c(19,"X")]
    drawsAB <- draws[,c(2:3,4,5,8)]
    drawsABC <- draws[,c(2:3,4,5,8,19,"X")]
    drawsBACA <- draws[,c(4,5,8,2,3,19,"X",2,3)]

    expect_equal(cbind(drawsA, drawsB), drawsAB)
    expect_equal(cbind(drawsA, drawsB, drawsC), drawsABC)
    expect_equal(cbind(drawsB, drawsA, drawsC, drawsA), drawsBACA)

})

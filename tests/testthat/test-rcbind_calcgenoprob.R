context("rbind and cbind calc_genoprob")

test_that("rbind.calc_genoprob works for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    probsA <- probs[1:5,]
    probsB <- probs[6:12,]
    probsC <- probs[13:20,]
    probsAB <- probs[1:12,]
    probsABC <- probs[1:20,]
    probsBACA <- probs[c(6:12, 1:5, 13:20, 1:5),]

    expect_equal(rbind(probsA, probsB), probsAB)
    expect_equal(rbind(probsA, probsB, probsC), probsABC)
    expect_equal(rbind(probsB, probsA, probsC, probsA), probsBACA)

})


test_that("rbind.calc_genoprob works for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    probsA <- probs[2:20,]
    probsB <- probs[41:60,]
    probsC <- probs[102:201,]
    probsAB <- probs[c(2:20,41:60),]
    probsABC <- probs[c(2:20,41:60,102:201),]
    probsBACA <- probs[c(41:60,2:20,102:201,2:20),]

    expect_equal(rbind(probsA, probsB), probsAB)
    expect_equal(rbind(probsA, probsB, probsC), probsABC)
    expect_equal(rbind(probsB, probsA, probsC, probsA), probsBACA)

})

test_that("cbind.calc_genoprob works for grav", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2[1:10,], map, error_prob=0.002)
    probsA <- probs[,1:2]
    probsB <- probs[,5]
    probsC <- probs[,3:4]
    probsAB <- probs[,c(1:2,5)]
    probsABC <- probs[,c(1,2,5,3,4)]
    probsBACA <- probs[,c(5,1:2,3:4,1:2)]

    expect_equal(cbind(probsA, probsB), probsAB)
    expect_equal(cbind(probsA, probsB, probsC), probsABC)
    expect_equal(cbind(probsB, probsA, probsC, probsA), probsBACA)

})


test_that("cbind.calc_genoprob works for iron", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron[6:21,], map, error_prob=0.002)
    probsA <- probs[,2:3]
    probsB <- probs[,c(4,5,8)]
    probsC <- probs[,c(19,"X")]
    probsAB <- probs[,c(2:3,4,5,8)]
    probsABC <- probs[,c(2:3,4,5,8,19,"X")]
    probsBACA <- probs[,c(4,5,8,2,3,19,"X",2,3)]

    expect_equal(cbind(probsA, probsB), probsAB)
    expect_equal(cbind(probsA, probsB, probsC), probsABC)
    expect_equal(cbind(probsB, probsA, probsC, probsA), probsBACA)

})

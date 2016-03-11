context("cbind and rbind for scan1 objects")

test_that("cbind.scan1() works for scan1() results", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

    out1 <- scan1(probs, grav2$pheno[,1,drop=FALSE])
    out2 <- scan1(probs, grav2$pheno[,2,drop=FALSE])
    out12 <- scan1(probs, grav2$pheno[,1:2])

    expect_equal(cbind(out1, out2), out12)

})

test_that("cbind.scan1() works for scan1/LMM results", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    k <- calc_kinship(probs)
    kloco <- calc_kinship(probs, "loco")

    out1 <- scan1(probs, grav2$pheno[,1,drop=FALSE], k)
    out2 <- scan1(probs, grav2$pheno[,2,drop=FALSE], k)
    out12 <- scan1(probs, grav2$pheno[,1:2], k)
    expect_equal(cbind(out1, out2), out12)

    out1 <- scan1(probs, grav2$pheno[,1,drop=FALSE], kloco)
    out2 <- scan1(probs, grav2$pheno[,2,drop=FALSE], kloco)
    out12 <- scan1(probs, grav2$pheno[,1:2], kloco)
    expect_equal(cbind(out1, out2), out12)

})

test_that("rbind.scan1() works for scan1() results", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

    phe <- grav2$pheno[,1,drop=FALSE]
    out1 <- scan1(probs[,1], phe)
    out2 <- scan1(probs[,2:3], phe)
    out3 <- scan1(probs[,5], phe)
    out12 <- scan1(probs[,1:3], phe)
    out123 <- scan1(probs[,c(1:3,5)], phe)
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe)

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

})


test_that("rbind.scan1() works for scan1() results with multiple columns", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

    phe <- grav2$pheno[,15:18,drop=FALSE]
    out1 <- scan1(probs[,1], phe)
    out2 <- scan1(probs[,2:3], phe)
    out3 <- scan1(probs[,5], phe)
    out12 <- scan1(probs[,1:3], phe)
    out123 <- scan1(probs[,c(1:3,5)], phe)
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe)

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

})


test_that("rbind.scan1() works for scan1() results", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    k <- calc_kinship(probs)
    kloco <- calc_kinship(probs, "loco")

    phe <- grav2$pheno[,1,drop=FALSE]
    out1 <- scan1(probs[,1], phe, k)
    out2 <- scan1(probs[,2:3], phe, k)
    out3 <- scan1(probs[,5], phe, k)
    out12 <- scan1(probs[,1:3], phe, k)
    out12$hsq <- rbind(out1$hsq, out2$hsq) # small adjustment
    out123 <- scan1(probs[,c(1:3,5)], phe, k)
    out123$hsq <- rbind(out1$hsq, out2$hsq, out3$hsq) # small adjustment
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe, k)
    out2123$hsq <- rbind(out2$hsq, out1$hsq, out2$hsq, out3$hsq) # small adjustment

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

    out1 <- scan1(probs[,1], phe, kloco[1])
    out2 <- scan1(probs[,2:3], phe, kloco[2:3])
    out3 <- scan1(probs[,5], phe, kloco[5])
    out12 <- scan1(probs[,1:3], phe, kloco[1:3])
    out123 <- scan1(probs[,c(1:3,5)], phe, kloco[c(1:3,5)])
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe, kloco[c(2:3,1,2:3,5)])

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

})


test_that("rbind.scan1() works for scan1() results with multiple columns", {

    library(qtl2geno)
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    k <- calc_kinship(probs)
    kloco <- calc_kinship(probs, "loco")

    phe <- grav2$pheno[,15:18,drop=FALSE]
    out1 <- scan1(probs[,1], phe, k)
    out2 <- scan1(probs[,2:3], phe, k)
    out3 <- scan1(probs[,5], phe, k)
    out12 <- scan1(probs[,1:3], phe, k)
    out12$hsq <- rbind(out1$hsq, out2$hsq) # small adjustment
    out123 <- scan1(probs[,c(1:3,5)], phe, k)
    out123$hsq <- rbind(out1$hsq, out2$hsq, out3$hsq) # small adjustment
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe, k)
    out2123$hsq <- rbind(out2$hsq, out1$hsq, out2$hsq, out3$hsq) # small adjustment

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

    out1 <- scan1(probs[,1], phe, kloco[1])
    out2 <- scan1(probs[,2:3], phe, kloco[2:3])
    out3 <- scan1(probs[,5], phe, kloco[5])
    out12 <- scan1(probs[,1:3], phe, kloco[1:3])
    out123 <- scan1(probs[,c(1:3,5)], phe, kloco[c(1:3,5)])
    out2123 <- scan1(probs[,c(2:3,1,2:3,5)], phe, kloco[c(2:3,1,2:3,5)])

    expect_equal(rbind(out1, out2), out12)
    expect_equal(rbind(out1, out2, out3), out123)
    expect_equal(rbind(out2, out1, out2, out3), out2123)

})

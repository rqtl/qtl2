context("calc_genoprob2")
suppressMessages(library(qtl))

test_that("backcross autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]

    hyper <- convert2cross2(hyper)
    pr <- calc_genoprob(hyper, error_prob=0.002, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, error_prob=0.002, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("intercross autosome", {

    data(listeria)
    chr <- c(1, 4, 14, 18)
    listeria <- listeria[chr,]

    listeria <- convert2cross2(listeria)
    map <- insert_pseudomarkers(listeria$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(listeria, map, error_prob=0.01, lowmem=TRUE)
    pr2 <- calc_genoprob(listeria, map, error_prob=0.01, lowmem=FALSE)

    expect_equal(pr, pr2)

})


test_that("risib autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.002, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.002, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("riself autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.002, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.002, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("f2 X chr", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(fake.f2, map, error_prob=0.01, lowmem=TRUE)
    pr2 <- calc_genoprob(fake.f2, map, error_prob=0.01, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("bc X chr", {

    data(hyper)
    hyper <- hyper["X",]
    # make it half female, half male
    hyper$pheno$sex <- rep(c("female", "male"), nind(hyper)/2)

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("f2 X chr all males", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 1

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(fake.f2, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(fake.f2, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})


test_that("f2 X chr all females", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(fake.f2, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(fake.f2, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("f2 X chr all females forw", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 0

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(fake.f2, map, error_prob=0.0001, lowmem=TRUE)
    pr2 <- calc_genoprob(fake.f2, map, error_prob=0.0001, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("f2 X chr all females rev", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 1

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(fake.f2, map, error_prob=0.0001, lowmem=TRUE)
    pr2 <- calc_genoprob(fake.f2, map, error_prob=0.0001, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("bc X chr all males", {

    data(hyper)
    hyper <- hyper["X",]

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})


test_that("bc X chr all females", {

    data(hyper)
    hyper <- hyper["X",]
    hyper$pheno$sex <- "female"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("doubled haploids", {

    data(hyper)
    hyper <- hyper[1,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})


test_that("haploids calc_genoprob", {

    data(hyper)
    hyper <- hyper[2,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=TRUE)
    pr2 <- calc_genoprob(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(pr, pr2)

})

test_that("backcross autosome with markers at same location", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    grav2 <- grav2[1:4,1]

    # put some markers at same location
    grav2$gmap[[1]][4] <- grav2$gmap[[1]][3]

    pr <- calc_genoprob(grav2, err=0, lowmem=TRUE)
    pr2 <- calc_genoprob(grav2, err=0, lowmem=FALSE)

    expect_true( all(!is.na(pr[[1]])) )
    expect_equal(pr, pr2)

})

test_that("calc_genoprob works when multi-core", {
    if(isnt_karl()) skip("this test only run locally")

    data(hyper)
    hyper2 <- convert2cross2(hyper)
    pr <- calc_genoprob(hyper2, error_prob=0.002, lowmem=TRUE)
    pr_mc <- calc_genoprob(hyper2, error_prob=0.002, cores=4, lowmem=TRUE)
    expect_equal(pr_mc, pr)

    pr2 <- calc_genoprob(hyper2, error_prob=0.002, lowmem=FALSE)
    expect_equal(pr2, pr)
    pr2_mc <- calc_genoprob(hyper2, error_prob=0.002, cores=4, lowmem=FALSE)
    expect_equal(pr2_mc, pr)

    data(listeria)
    listeria2 <- convert2cross2(listeria)
    map <- insert_pseudomarkers(listeria2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(listeria2, map, error_prob=0.01, lowmem=TRUE)
    pr_mc <- calc_genoprob(listeria2, map, error_prob=0.01, cores=4, lowmem=TRUE)
    expect_equal(pr_mc, pr)

    pr2 <- calc_genoprob(listeria2, map, error_prob=0.01, lowmem=FALSE)
    expect_equal(pr2, pr)
    pr2_mc <- calc_genoprob(listeria2, map, error_prob=0.01, cores=4, lowmem=FALSE)
    expect_equal(pr2_mc, pr)

})

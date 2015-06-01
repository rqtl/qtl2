
context("est_map")
library(qtl)

test_that("est_map for backcross autosome matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    newmap <- est.map(hyper, err=0.002, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    hyper2 <- convert2cross2(hyper)
    newmap2 <- est_map(hyper2, err=0.002, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

test_that("est_map for intercross autosome matches R/qtl", {

    data(listeria)
    chr <- c(4, 14, 18)
    listeria <- listeria[chr,]
    newmap <- est.map(listeria, err=0.01, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    listeria2 <- convert2cross2(listeria)
    newmap2 <- est_map(listeria2, err=0.01, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))


})

test_that("f2 X chr est_map matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    newmap <- est.map(fake.f2, err=0.01, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    fake.f2.2 <- convert2cross2(fake.f2)
    newmap2 <- est_map(fake.f2.2, err=0.01, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

test_that("bc X chr calc_genoprob matches R/qtl", {

    set.seed(19115167)
    xmap <- sim.map(100, n.mar=11, anchor.tel=TRUE, include.x=TRUE, eq.spacing=TRUE)
    n.ind <- 100
    cross <- sim.cross(xmap, n.ind=n.ind, type="bc",
                       error.prob=0.01, missing.prob=0.05)
    cross$pheno$sex <- rep(c(0,1), n.ind/2)
    sexpgm <- getsex(cross)

    newmap <- est.map(cross, err=0.01, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    cross2 <- convert2cross2(cross)
    newmap2 <- est_map(cross2, err=0.01, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

test_that("est_map for RIself matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    hyper2 <- convert2cross2(hyper)
    newmap2 <- est_map(hyper2, err=0.002, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})


test_that("est_map for RIsib matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    hyper2 <- convert2cross2(hyper)
    newmap2 <- est_map(hyper2, err=0.002, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

test_that("est_map for doubled haploids matches R/qtl", {

    data(hyper)
    hyper <- hyper[3,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"

    newmap <- est.map(hyper, err=0.002, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    hyper2 <- convert2cross2(hyper)
    newmap2 <- est_map(hyper2, err=0.002, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

test_that("est_map for haploids matches R/qtl", {

    data(hyper)
    hyper <- hyper[4,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"

    newmap <- est.map(hyper, err=0.002, tol=1e-8)
    newmap <- lapply(newmap, unclass)

    hyper2 <- convert2cross2(hyper)
    newmap2 <- est_map(hyper2, err=0.002, tol=1e-8)

    expect_equivalent(newmap, newmap2)
    expect_equal(lapply(newmap, attr, "loglik"),
                 lapply(newmap2, attr, "loglik"))

})

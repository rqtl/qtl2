
context("est_map")
library(qtl)

test_that("est_map for backcross autosome matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    newmap <- est.map(hyper, err=0.002, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:nchr(hyper)) {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("bc", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})

test_that("est_map for intercross autosome matches R/qtl", {

    data(listeria)
    chr <- c(4, 14, 18)
    listeria <- listeria[chr,]
    newmap <- est.map(listeria, err=0.01, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:nchr(listeria)) {
        g <- listeria$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- listeria$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("f2", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.01, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})

test_that("f2 X chr est_map matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]

    # omit individuals with <2 genotypes
    fake.f2 <- fake.f2[,ntyped(fake.f2) > 1]

    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    gorig <- g

    # males -> 1/3
    gmale <- g[sexpgm$sex==1,]
    gmale[gmale==2] <- 3
    g[sexpgm$sex==1,] <- gmale

    # female rev -> 2/3
    gfemalerev <- g[sexpgm$sex==0 & sexpgm$pgm==1,]
    gfemalerev[gfemalerev==1] <- 3
    g[sexpgm$sex==0 & sexpgm$pgm==1,] <- gfemalerev

    # R/qtl est.map
    newmap <- est.map(fake.f2, err=0.01, tol=1e-8)
    rf_qtl <- as.numeric(mf.h(diff(newmap[[1]])))

    map <- fake.f2$geno[["X"]]$map
    rf <- mf.h(diff(map))

    rf_qtl2 <- est_map("f2", t(g), TRUE, (sexpgm$sex==0), rbind(sexpgm$pgm),
                       rf, 0.01, 10000, 1e-8, FALSE)

    # treated as bc autosome
    rf_qtl2_alt <- est_map("bc", t(gorig), FALSE, rep(TRUE, nind(fake.f2)), rbind(sexpgm$pgm),
                           rf, 0.01, 10000, 1e-8, FALSE)

    expect_equal(rf_qtl2, rf_qtl2_alt)
    expect_equal(attr(rf_qtl2, "loglik"), attr(rf_qtl2_alt, "loglik"))

    # same as R/qtl?
    expect_equivalent(rf_qtl, rf_qtl2)
    expect_equal(attr(newmap[[1]], "loglik"), attr(rf_qtl2, "loglik"))
})

test_that("bc X chr calc_genoprob matches R/qtl", {

    set.seed(19115167)
    xmap <- sim.map(100, n.mar=11, anchor.tel=TRUE, include.x=TRUE, eq.spacing=TRUE)
    n.ind <- 100
    cross <- sim.cross(xmap, n.ind=n.ind, type="bc",
                       error.prob=0.01, missing.prob=0.05)
    cross$pheno$sex <- rep(c(0,1), n.ind/2)
    sexpgm <- getsex(cross)

    g <- cross$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    gorig <- g

    # males -> 1/3
    gmale <- g[sexpgm$sex==1,]
    gmale[gmale==2] <- 3
    g[sexpgm$sex==1,] <- gmale

    # R/qtl est.map
    newmap <- est.map(cross, err=0.01, tol=1e-8)
    rf_qtl <- as.numeric(mf.h(diff(newmap[[1]])))

    map <- cross$geno[["X"]]$map
    rf <- mf.h(diff(map))

    rf_qtl2 <- est_map("bc", t(g), TRUE, (sexpgm$sex==0),
                       matrix(nrow=0, ncol=nind(cross)),
                       rf, 0.01, 10000, 1e-8, FALSE)

    # treat as autosome
    rf_qtl2_alt <- est_map("bc", t(gorig), FALSE, (sexpgm$sex==0),
                           matrix(nrow=0, ncol=nind(cross)),
                           rf, 0.01, 10000, 1e-8, FALSE)

    expect_equal(rf_qtl2, rf_qtl2_alt)
    expect_equal(attr(rf_qtl2, "loglik"), attr(rf_qtl2_alt, "loglik"))

    # same as R/qtl?
    expect_equivalent(rf_qtl, rf_qtl2)
    expect_equal(attr(newmap[[1]], "loglik"), attr(rf_qtl2, "loglik"))

})

test_that("est_map for RIself matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:nchr(hyper)) {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("riself", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})


test_that("est_map for RIsib matches R/qtl", {

    data(hyper)
    chr <- c(3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:nchr(hyper)) {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("risib", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})

test_that("est_map for doubled haploids matches R/qtl", {

    data(hyper)
    hyper <- hyper[3,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in "3") {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("dh", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})

test_that("est_map for haploids matches R/qtl", {

    data(hyper)
    hyper <- hyper[4,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"
    newmap <- est.map(hyper, err=0.002, tol=1e-8)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in "4") {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("haploid", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-8, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})


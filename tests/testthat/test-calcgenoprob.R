
context("calc_genoprob")

test_that("backcross autosome calc_genoprob matches R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[1:19,]
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)

    for(i in 1:19) {
       g <- hyper$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- hyper$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("bc", t(g), FALSE, rep(FALSE, nind(hyper)),
                                matrix(ncol=nind(hyper), nrow=0),
                                rf, map_index, 0.002)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})

test_that("intercross autosome calc_genoprob matches R/qtl", {

    library(qtl)
    data(listeria)
    listeria <- listeria[1:19,]
    listeria <- calc.genoprob(listeria, step=1, stepwidth="max", err=0.01)

    for(i in 1:19) {
       g <- listeria$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- listeria$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("f2", t(g), FALSE, rep(FALSE, nind(listeria)),
                                matrix(ncol=nind(listeria), nrow=0),
                                rf, map_index, 0.01)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})


test_that("risib autosome calc_genoprob matches R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[1:19,]
    class(hyper)[1] <- "risib"
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)

    for(i in 1:19) {
       g <- hyper$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- hyper$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("risib", t(g), FALSE, rep(FALSE, nind(hyper)),
                                matrix(ncol=nind(hyper), nrow=0),
                                rf, map_index, 0.002)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})

test_that("riself autosome calc_genoprob matches R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[1:19,]
    class(hyper)[1] <- "riself"
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)

    for(i in 1:19) {
       g <- hyper$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- hyper$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("riself", t(g), FALSE, rep(FALSE, nind(hyper)),
                                matrix(ncol=nind(hyper), nrow=0),
                                rf, map_index, 0.002)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})

test_that("f2 X chr calc_genoprob matches R/qtl", {

    library(qtl)
    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    # males -> 1/3
    gmale <- g[sexpgm$sex==1,]
    gmale[gmale==2] <- 3
    g[sexpgm$sex==1,] <- gmale

    # female rev -> 2/3
    gfemalerev <- g[sexpgm$sex==0 & sexpgm$pgm==1,]
    gfemalerev[gfemalerev==1] <- 3
    g[sexpgm$sex==0 & sexpgm$pgm==1,] <- gfemalerev

    # R/qtl calc.genoprob
    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.01)

    pr <- fake.f2$geno[["X"]]$prob
    map <- attr(pr, "map")
    rf <- mf.h(diff(map))
    map_index <- match(names(map), colnames(g))-1
    map_index[is.na(map_index)] <- -1

    pr_qtl2 <- calc_genoprob("f2", t(g), TRUE, (sexpgm$sex==0), rbind(sexpgm$pgm),
                             rf, map_index, 0.01)
    pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

    pr_qtl <- reviseXdata("f2", "simple", sexpgm=sexpgm, prob=pr,
                          cross.attr=attributes(fake.f2))

    expect_equivalent(pr_qtl, pr_qtl2)

})

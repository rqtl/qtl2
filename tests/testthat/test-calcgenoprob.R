
context("calc_genoprob")
library(qtl)

test_that("backcross autosome calc_genoprob matches R/qtl", {

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

    pr_qtl <- reviseXdata("f2", "full", sexpgm=sexpgm, prob=pr,
                          cross.attr=attributes(fake.f2))

    expect_equivalent(pr_qtl, pr_qtl2)

})

test_that("bc X chr calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    # make it half female, half male
    hyper$pheno$sex <- rep(c("female", "male"), nind(hyper)/2)
    sexpgm <- getsex(hyper)

    g <- hyper$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    # males -> 1/3
    gmale <- g[sexpgm$sex==1,]
    gmale[gmale==2] <- 3
    g[sexpgm$sex==1,] <- gmale

    # R/qtl calc.genoprob
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)

    pr <- hyper$geno[["X"]]$prob
    map <- attr(pr, "map")
    rf <- mf.h(diff(map))
    map_index <- match(names(map), colnames(g))-1
    map_index[is.na(map_index)] <- -1

    pr_qtl2 <- calc_genoprob("bc", t(g), TRUE, (sexpgm$sex==0),
                             matrix(nrow=0, ncol=nind(hyper)),
                             rf, map_index, 0.02)
    pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

    pr_qtl <- reviseXdata("bc", "full", sexpgm=sexpgm, prob=pr,
                          cross.attr=attributes(hyper))

    expect_equivalent(pr_qtl, pr_qtl2)

})

test_that("f2 X chr all males calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 1
    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    # males -> 1/3
    g[g==2] <- 3

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

    # initial four probs all == 0
    for(i in 1:4)
        expect_true(all(pr_qtl2[,,i] == 0))

    expect_equivalent(pr, pr_qtl2[,,5:6])

})


test_that("f2 X chr all females calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

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

    # last two probs all 0
    for(i in 5:6)
        expect_true(all(pr_qtl2[,,i]==0))

    pr_qtl <- reviseXdata("f2", "full", sexpgm=sexpgm, prob=pr,
                          cross.attr=attributes(fake.f2))

    expect_equivalent(pr_qtl, pr_qtl2[,,1:4])

})

test_that("f2 X chr all females forw calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 0
    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

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

    # expect last four probs == 0
    for(i in 3:6)
        expect_true(all(pr_qtl2[,,i] == 0))

    expect_equivalent(pr, pr_qtl2[,,1:2])

})

test_that("f2 X chr all females rev calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 1
    sexpgm <- getsex(fake.f2)

    g <- fake.f2$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

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

    # first two and last two probs all 0
    for(i in c(1:2, 5:6))
        expect_true(all(pr_qtl2[,,i]==0))

    expect_equivalent(pr, pr_qtl2[,,4:3])

})

test_that("bc X chr all males calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    sexpgm <- getsex(hyper)

    g <- hyper$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    # males -> 1/3
    g[g==2] <- 3

    # R/qtl calc.genoprob
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)

    pr <- hyper$geno[["X"]]$prob
    map <- attr(pr, "map")
    rf <- mf.h(diff(map))
    map_index <- match(names(map), colnames(g))-1
    map_index[is.na(map_index)] <- -1

    pr_qtl2 <- calc_genoprob("bc", t(g), TRUE, (sexpgm$sex==0),
                             matrix(nrow=0, ncol=nind(hyper)),
                             rf, map_index, 0.02)
    pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

    # first two probs all 0
    for(i in 1:2)
        expect_true(all(pr_qtl2[,,i]==0))

    expect_equivalent(pr, pr_qtl2[,,3:4])

})


test_that("bc X chr all females calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    hyper$pheno$sex <- "female"
    sexpgm <- getsex(hyper)

    g <- hyper$geno[["X"]]$data
    g[is.na(g)] <- 0 # missing -> 0

    # R/qtl calc.genoprob
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)

    pr <- hyper$geno[["X"]]$prob
    map <- attr(pr, "map")
    rf <- mf.h(diff(map))
    map_index <- match(names(map), colnames(g))-1
    map_index[is.na(map_index)] <- -1

    pr_qtl2 <- calc_genoprob("bc", t(g), TRUE, (sexpgm$sex==0),
                             matrix(nrow=0, ncol=nind(hyper)),
                             rf, map_index, 0.02)
    pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

    # last two probs all 0
    for(i in 3:4)
        expect_true(all(pr_qtl2[,,i]==0))

    expect_equivalent(pr, pr_qtl2[,,1:2])

})

test_that("doubled haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[1,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)

    for(i in "1") {
       g <- hyper$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- hyper$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("dh", t(g), FALSE, rep(FALSE, nind(hyper)),
                                matrix(ncol=nind(hyper), nrow=0),
                                rf, map_index, 0.002)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})


test_that("haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[2,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"
    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)

    for(i in "2") {
       g <- hyper$geno[[i]]$data
       g[is.na(g)] <- 0
       pr <- hyper$geno[[i]]$prob
       map <- attr(pr, "map")
       rf <- mf.h(diff(map))
       map_index <- match(names(map), colnames(g))-1
       map_index[is.na(map_index)] <- -1

       pr_qtl2 <- calc_genoprob("haploid", t(g), FALSE, rep(FALSE, nind(hyper)),
                                matrix(ncol=nind(hyper), nrow=0),
                                rf, map_index, 0.002)
       pr_qtl2 <- aperm(pr_qtl2, c(2,3,1))

       expect_equivalent(pr, pr_qtl2)
   }

})



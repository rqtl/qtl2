context("calc_genoprob")
suppressMessages(library(qtl))

grab_prob_rqtl <- function(cross) lapply(cross$geno, function(a) aperm(a$prob, c(1,3,2)))


test_that("backcross autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]

    hyper <- calc.genoprob(hyper, err=0.002)
    pr <- grab_prob_rqtl(hyper)

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, error_prob=0.002)

    expect_equivalent(pr, pr2, tol=2e-5)

})

test_that("intercross autosome calc_genoprob matches R/qtl", {

    data(listeria)
    chr <- c(1, 4, 14, 18)
    listeria <- listeria[chr,]

    listeria <- calc.genoprob(listeria, step=1, stepwidth="max", err=0.01)
    pr <- grab_prob_rqtl(listeria)

    listeria2 <- convert2cross2(listeria)
    map <- insert_pseudomarkers(listeria2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(listeria2, map, error_prob=0.01)

    expect_equivalent(pr, pr2)

})


test_that("risib autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)
    pr <- grab_prob_rqtl(hyper)

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.002)

    expect_equivalent(pr, pr2, tol=1e-5)

})

test_that("riself autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)
    pr <- grab_prob_rqtl(hyper)

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.002)

    expect_equivalent(pr, pr2, tol=1e-5)

})

test_that("f2 X chr calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.01)
    pr <- fake.f2$geno[[1]]$prob
    pr <- reviseXdata("f2", "full", getsex(fake.f2), prob=pr, cross.attr=attributes(fake.f2))
    pr <- aperm(pr, c(1,3,2))

    fake.f2.2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2.2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(fake.f2.2, map, error_prob=0.01)[[1]]

    expect_equivalent(pr, pr2)

})

test_that("bc X chr calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    # make it half female, half male
    hyper$pheno$sex <- rep(c("female", "male"), nind(hyper)/2)

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[["X"]]$prob
    pr <- reviseXdata("bc", "full", getsex(hyper), prob=pr, cross.attr=attributes(hyper))
    pr <- aperm(pr, c(1,3,2))

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all males calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 1

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.02)
    pr <- fake.f2$geno[["X"]]$prob
    pr <- aperm(pr, c(1,3,2))

    fake.f2.2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2.2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(fake.f2.2, map, error_prob=0.02)[[1]][,5:6,]

    expect_equivalent(pr, pr2)

})


test_that("f2 X chr all females calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.02)
    pr <- fake.f2$geno[["X"]]$prob
    pr <- aperm(pr, c(1,3,2))

    fake.f2.2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2.2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(fake.f2.2, map, error_prob=0.02)[[1]]
    pr2 <- pr2[,1:2,] + pr2[,4:3,] # recode as in R/qtl

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all females forw calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 0

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.0001)
    pr <- fake.f2$geno[["X"]]$prob
    pr <- aperm(pr, c(1,3,2))

    fake.f2.2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2.2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(fake.f2.2, map, error_prob=0.0001)[[1]][,1:2,]

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all females rev calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 1
    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.0001)
    pr <- fake.f2$geno[["X"]]$prob
    pr <- aperm(pr, c(1,3,2))

    fake.f2.2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2.2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(fake.f2.2, map, error_prob=0.0001)[[1]][,4:3,]

    expect_equivalent(pr, pr2)

})

test_that("bc X chr all males calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    sexpgm <- getsex(hyper)

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[["X"]]$prob
    pr <- reviseXdata("bc", "full", getsex(hyper), prob=pr, cross.attr=attributes(hyper))
    pr <- aperm(pr, c(1,3,2))

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.02)[[1]][,3:4,]

    expect_equivalent(pr, pr2)

})


test_that("bc X chr all females calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    hyper$pheno$sex <- "female"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[["X"]]$prob
    pr <- reviseXdata("bc", "full", getsex(hyper), prob=pr, cross.attr=attributes(hyper))
    pr <- aperm(pr, c(1,3,2))

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.02)[[1]]
    pr2 <- pr2[,1:2,] + pr2[,3:4,]

    expect_equivalent(pr, pr2)

})

test_that("doubled haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[1,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[[1]]$prob
    pr <- aperm(pr, c(1,3,2))

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2, tol=1e-6)

})


test_that("haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[2,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[[1]]$prob
    pr <- aperm(pr, c(1,3,2))

    hyper2 <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper2$gmap, step=1, stepwidth="max")
    pr2 <- calc_genoprob(hyper2, map, error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2)

})

test_that("backcross autosome calc_genoprob with markers at same location", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[1:4,1]

    # put some markers at same location
    grav2$gmap[[1]][4] <- grav2$gmap[[1]][3]

    pr <- calc_genoprob(grav2, err=0)

    expect_true( all(!is.na(pr[[1]])) )

})

test_that("calc_genoprob works when multi-core", {
    if(isnt_karl()) skip("this test only run locally")

    data(hyper)
    hyper2 <- convert2cross2(hyper)
    pr <- calc_genoprob(hyper2, error_prob=0.002)
    pr_mc <- calc_genoprob(hyper2, error_prob=0.002, cores=4)
    expect_equal(pr_mc, pr)


    data(listeria)
    listeria2 <- convert2cross2(listeria)
    map <- insert_pseudomarkers(listeria2$gmap, step=1, stepwidth="max")
    pr <- calc_genoprob(listeria2, map, error_prob=0.01)
    pr_mc <- calc_genoprob(listeria2, map, error_prob=0.01, cores=4)

    expect_equal(pr_mc, pr)

})

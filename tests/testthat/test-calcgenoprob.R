
context("calc_genoprob")
library(qtl)

test_that("backcross autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]

    hyper <- calc.genoprob(hyper, err=0.002)
    pr <- lapply(hyper$geno, "[[", "prob")

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, error_prob=0.002)

    expect_equivalent(pr, pr2)

})

test_that("intercross autosome calc_genoprob matches R/qtl", {

    data(listeria)
    chr <- c(1, 4, 14, 18)
    listeria <- listeria[chr,]

    listeria <- calc.genoprob(listeria, step=1, stepwidth="max", err=0.01)
    pr <- lapply(listeria$geno, "[[", "prob")

    listeria2 <- convert2cross2(listeria)
    pr2 <- calc_genoprob(listeria2, step=1, stepwidth="max", error_prob=0.01)

    expect_equivalent(pr, pr2)

})


test_that("risib autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)
    pr <- lapply(hyper$geno, "[[", "prob")

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.002)

    expect_equivalent(pr, pr2)

})

test_that("riself autosome calc_genoprob matches R/qtl", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", err=0.002)
    pr <- lapply(hyper$geno, "[[", "prob")

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.002)

    expect_equivalent(pr, pr2)


})

test_that("f2 X chr calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.01)
    pr <- fake.f2$geno[[1]]$prob
    pr <- reviseXdata("f2", "full", getsex(fake.f2), prob=pr, cross.attr=attributes(fake.f2))

    fake.f2.2 <- convert2cross2(fake.f2)
    pr2 <- calc_genoprob(fake.f2.2, step=1, stepwidth="max", error_prob=0.01)[[1]]

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

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all males calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 1

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.02)
    pr <- fake.f2$geno[["X"]]$prob

    fake.f2.2 <- convert2cross2(fake.f2)
    pr2 <- calc_genoprob(fake.f2.2, step=1, stepwidth="max", error_prob=0.02)[[1]][,,5:6]

    expect_equivalent(pr, pr2)

})


test_that("f2 X chr all females calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.02)
    pr <- fake.f2$geno[["X"]]$prob

    fake.f2.2 <- convert2cross2(fake.f2)
    pr2 <- calc_genoprob(fake.f2.2, step=1, stepwidth="max", error_prob=0.02)[[1]]
    pr2 <- pr2[,,1:2] + pr2[,,4:3] # recode as in R/qtl

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all females forw calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 0

    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.0001)
    pr <- fake.f2$geno[["X"]]$prob

    fake.f2.2 <- convert2cross2(fake.f2)
    pr2 <- calc_genoprob(fake.f2.2, step=1, stepwidth="max", error_prob=0.0001)[[1]][,,1:2]

    expect_equivalent(pr, pr2)

})

test_that("f2 X chr all females rev calc_genoprob matches R/qtl", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 1
    fake.f2 <- calc.genoprob(fake.f2, step=1, stepwidth="max", error.prob=0.0001)
    pr <- fake.f2$geno[["X"]]$prob

    fake.f2.2 <- convert2cross2(fake.f2)
    pr2 <- calc_genoprob(fake.f2.2, step=1, stepwidth="max", error_prob=0.0001)[[1]][,,4:3]

    expect_equivalent(pr, pr2)

})

test_that("bc X chr all males calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    sexpgm <- getsex(hyper)

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[["X"]]$prob
    pr <- reviseXdata("bc", "full", getsex(hyper), prob=pr, cross.attr=attributes(hyper))

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.02)[[1]][,,3:4]

    expect_equivalent(pr, pr2)

})


test_that("bc X chr all females calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper["X",]
    hyper$pheno$sex <- "female"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[["X"]]$prob
    pr <- reviseXdata("bc", "full", getsex(hyper), prob=pr, cross.attr=attributes(hyper))

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.02)[[1]]
    pr2 <- pr2[,,1:2] + pr2[,,3:4]

    expect_equivalent(pr, pr2)

})

test_that("doubled haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[1,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[[1]]$prob

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2)

})


test_that("haploids calc_genoprob matches R/qtl", {

    data(hyper)
    hyper <- hyper[2,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"

    hyper <- calc.genoprob(hyper, step=1, stepwidth="max", error.prob=0.02)
    pr <- hyper$geno[[1]]$prob

    hyper2 <- convert2cross2(hyper)
    pr2 <- calc_genoprob(hyper2, step=1, stepwidth="max", error_prob=0.02)[[1]]

    expect_equivalent(pr, pr2)

})



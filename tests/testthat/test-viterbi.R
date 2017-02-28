context("viterbi")
suppressMessages(library(qtl))

test_that("backcross autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]

    hyper <- convert2cross2(hyper)
    set.seed(20150524)
    g <- viterbi(hyper, error_prob=0.002, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, error_prob=0.002, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("intercross autosome", {

    data(listeria)
    chr <- c(1, 4, 14, 18)
    listeria <- listeria[chr,]

    listeria <- convert2cross2(listeria)
    map <- insert_pseudomarkers(listeria$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(listeria, map, error_prob=0.01, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(listeria, map, error_prob=0.01, lowmem=FALSE)

    expect_equal(g, g2)

})


test_that("risib autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "risib"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.002, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.002, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("riself autosome", {

    data(hyper)
    chr <- c(1, 3, 4, 17, 19)
    hyper <- hyper[chr,]
    class(hyper)[1] <- "riself"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.002, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.002, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("f2 X chr", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]

    fake.f2 <- convert2cross2(fake.f2)

    # order individuals by sex/cross
    sexcr <- paste(fake.f2$is_female, apply(fake.f2$cross_info, 1, paste, collapse=":"), sep=":")
    o <- order(sexcr)
    fake.f2 <- fake.f2[o,]

    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")

    set.seed(20150524)
    g <- viterbi(fake.f2, map, error_prob=0.01, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(fake.f2, map, error_prob=0.01, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("bc X chr", {

    data(hyper)
    hyper <- hyper["X",]
    # make it half female, half male
    hyper$pheno$sex <- rep(c("female", "male"), nind(hyper)/2)

    hyper <- convert2cross2(hyper)

    # order individuals by sex/cross
    sexcr <- paste(hyper$is_female, apply(hyper$cross_info, 1, paste, collapse=":"), sep=":")
    o <- order(sexcr)
    hyper <- hyper[o,]

    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")

    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("f2 X chr all males", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 1

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(fake.f2, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(fake.f2, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})


test_that("f2 X chr all females", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(fake.f2, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(fake.f2, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("f2 X chr all females forw", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 0

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(fake.f2, map, error_prob=0.0001, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(fake.f2, map, error_prob=0.0001, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("f2 X chr all females rev", {

    data(fake.f2)
    fake.f2 <- fake.f2["X",]
    fake.f2$pheno$sex <- 0
    fake.f2$pheno$pgm <- 1

    fake.f2 <- convert2cross2(fake.f2)
    map <- insert_pseudomarkers(fake.f2$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(fake.f2, map, error_prob=0.0001, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(fake.f2, map, error_prob=0.0001, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("bc X chr all males", {

    data(hyper)
    hyper <- hyper["X",]

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})


test_that("bc X chr all females", {

    data(hyper)
    hyper <- hyper["X",]
    hyper$pheno$sex <- "female"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("doubled haploids", {

    data(hyper)
    hyper <- hyper[1,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "dh"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})


test_that("haploids viterbi", {

    data(hyper)
    hyper <- hyper[2,]
    hyper$pheno <- hyper$pheno[,1,drop=FALSE]
    class(hyper)[1] <- "haploid"

    hyper <- convert2cross2(hyper)
    map <- insert_pseudomarkers(hyper$gmap, step=1, stepwidth="max")
    set.seed(20150524)
    g <- viterbi(hyper, map, error_prob=0.02, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(hyper, map, error_prob=0.02, lowmem=FALSE)

    expect_equal(g, g2)

})

test_that("backcross autosome with markers at same location", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    grav2 <- grav2[1:4,1]

    # put some markers at same location
    grav2$gmap[[1]][4] <- grav2$gmap[[1]][3]

    set.seed(20150524)
    g <- viterbi(grav2, error_prob=0, lowmem=TRUE)
    set.seed(20150524)
    g2 <- viterbi(grav2, error_prob=0, lowmem=FALSE)

    expect_true( all(!is.na(g[[1]])) )
    expect_equal(g, g2)

})

test_that("R/qtl2 gives same viterbi results as R/qtl when step=0, for backcross", {

    library(qtl)
    data(hyper)

    # step=0, error_prob=0, chromosomes where everyone has *some* data
    for(chr in c(1, 4, 7, 17)) {
        hypersub <- hyper[chr,]
        hyper2sub <- convert2cross2(hypersub)
        hypersub <- argmax.geno(hypersub, step=0, error.prob=0)
        agm <- hypersub$geno[[1]]$argmax
        agm2 <- viterbi(hyper2sub, error_prob=0)[[1]]
        expect_equivalent(agm, agm2)

        # lowmem version
        agm2b <- viterbi(hyper2sub, error_prob=0, lowmem=TRUE)[[1]]
        expect_equivalent(agm, agm2b)
    }

    # step=0, error_prob=0.002, chromosomes where everyone has *some* data
    for(chr in c(1, 4, 7, 17)) {
        hypersub <- hyper[chr,]
        hyper2sub <- convert2cross2(hypersub)
        hypersub <- argmax.geno(hypersub, step=0, error.prob=0.002)
        agm <- hypersub$geno[[1]]$argmax
        agm2 <- viterbi(hyper2sub, error_prob=0.002)[[1]]
        expect_equivalent(agm, agm2)

        # lowmem version
        agm2b <- viterbi(hyper2sub, error_prob=0.002, lowmem=TRUE)[[1]]
        expect_equivalent(agm, agm2b)
    }

})


test_that("R/qtl2 gives same viterbi results as R/qtl when step=0, for intercross", {

    library(qtl)
    data(listeria)

    # step=0, error_prob=0, chromosomes where everyone has *some* data
    for(chr in c(1, 4, 7, 17)) {
        listeriasub <- listeria[chr,]
        listeria2sub <- convert2cross2(listeriasub)
        listeriasub <- argmax.geno(listeriasub, step=0, error.prob=0)
        agm <- listeriasub$geno[[1]]$argmax
        agm2 <- viterbi(listeria2sub, error_prob=0)[[1]]
        expect_equivalent(agm, agm2)

        # lowmem version
        agm2b <- viterbi(listeria2sub, error_prob=0, lowmem=TRUE)[[1]]
        expect_equivalent(agm, agm2b)
    }

    # step=0, error_prob=0.002, chromosomes where everyone has *some* data
    for(chr in c(1, 4, 7, 17)) {
        listeriasub <- listeria[chr,]
        listeria2sub <- convert2cross2(listeriasub)
        listeriasub <- argmax.geno(listeriasub, step=0, error.prob=0.002)
        agm <- listeriasub$geno[[1]]$argmax
        agm2 <- viterbi(listeria2sub, error_prob=0.002)[[1]]
        expect_equivalent(agm, agm2)

        # lowmem version
        agm2b <- viterbi(listeria2sub, error_prob=0.002, lowmem=TRUE)[[1]]
        expect_equivalent(agm, agm2b)
    }

})

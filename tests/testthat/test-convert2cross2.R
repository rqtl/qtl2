
context("convert2cross2")
library(qtl)

test_that("convert2cross2 works appropriately for hyper data", {

    data(hyper)
    hyper2 <- convert2cross2(hyper)

    ids <- as.character(1:nind(hyper))

    # class and crosstype
    expect_equal(class(hyper2), "cross2")
    expect_equal(hyper2$crosstype, "bc")

    # X chr
    is_x <- rep(FALSE, nchr(hyper))
    names(is_x) <- c(1:19, "X")
    is_x["X"] <- TRUE
    expect_equal(hyper2$is_x_chr, is_x)

    # gmap
    gmap <- lapply(pull.map(hyper), function(a) { class(a) <- "numeric"; a })
    expect_equal(hyper2$gmap, gmap)

    # sex
    sex <- rep(1, nind(hyper))
    names(sex) <- ids
    expect_equal(hyper2$sex, sex)

    # cross_info
    cross_info <- matrix(ncol=0, nrow=nind(hyper))
    rownames(cross_info) <- ids
    expect_equal(hyper2$cross_info, cross_info)

    # pheno
    phe <- hyper$pheno[,1,drop=FALSE]
    rownames(phe) <- ids
    phe <- as.matrix(phe)
    expect_equal(hyper2$pheno, phe)

    # covar
    covar <- hyper$pheno[,2,drop=FALSE]
    rownames(covar) <- ids
    expect_equal(hyper2$covar, covar)

})

test_that("convert2cross2 works appropriately for fake.f2 data", {

    data(fake.f2)
    fake.f2.2 <- convert2cross2(fake.f2)

    ids <- as.character(1:nind(fake.f2))

    # class and crosstype
    expect_equal(class(fake.f2.2), "cross2")
    expect_equal(fake.f2.2$crosstype, "f2")

    # X chr
    is_x <- rep(FALSE, nchr(fake.f2))
    names(is_x) <- c(1:19, "X")
    is_x["X"] <- TRUE
    expect_equal(fake.f2.2$is_x_chr, is_x)

    # gmap
    gmap <- lapply(pull.map(fake.f2), function(a) { class(a) <- "numeric"; a })
    expect_equal(fake.f2.2$gmap, gmap)

    # sex
    sexpgm <- getsex(fake.f2)
    sex <- sexpgm$sex
    names(sex) <- ids
    expect_equal(fake.f2.2$sex, sex)

    # cross_info
    cross_info <- as.matrix(sexpgm$pgm)
    rownames(cross_info) <- ids
    expect_equal(fake.f2.2$cross_info, cross_info)

    # pheno
    phe <- fake.f2$pheno[,1,drop=FALSE]
    rownames(phe) <- ids
    phe <- as.matrix(phe)
    expect_equal(fake.f2.2$pheno, phe)

    # covar
    covar <- fake.f2$pheno[,2:3,drop=FALSE]
    rownames(covar) <- ids
    expect_equal(fake.f2.2$covar, covar)

})


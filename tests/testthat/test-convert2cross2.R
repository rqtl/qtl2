
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

    # genotypes
    for(i in seq(along=hyper$geno)) {
        go <- hyper$geno[[i]]$data
        if(is_x[i])
            go <- reviseXdata("bc", "simple", getsex(hyper),
                              geno=go, cross.attr=attributes(hyper),
                              force=TRUE)
        gn <- hyper2$geno[[i]]
        expect_true(all(is.na(go) | go > 0))
        go[is.na(go)] <- 0
        rownames(go) <- ids
        expect_equal(gn, go)
    }

    # gmap
    gmap <- lapply(pull.map(hyper),
                   function(a) { class(a) <- "numeric"; a })
    expect_equal(hyper2$gmap, gmap)

    # sex
    is_female <- rep(FALSE, nind(hyper))
    names(is_female) <- ids
    expect_equal(hyper2$is_female, is_female)

    # cross_info
    cross_info <- matrix(0L, ncol=0, nrow=nind(hyper))
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

    # genotypes
    for(i in seq(along=fake.f2.2$geno)) {
        go <- fake.f2$geno[[i]]$data
        if(is_x[i])
            go <- reviseXdata("f2", "simple", getsex(fake.f2),
                              geno=go, cross.attr=attributes(fake.f2))
        gn <- fake.f2.2$geno[[i]]
        expect_true(all(is.na(go) | go > 0))
        go[is.na(go)] <- 0
        rownames(go) <- ids
        expect_equal(gn, go)
    }

    # gmap
    gmap <- lapply(pull.map(fake.f2), function(a) { class(a) <- "numeric"; a })
    expect_equal(fake.f2.2$gmap, gmap)

    # sex
    sexpgm <- getsex(fake.f2)
    is_female <- (sexpgm$sex == 0)
    names(is_female) <- ids
    expect_equal(fake.f2.2$is_female, is_female)

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

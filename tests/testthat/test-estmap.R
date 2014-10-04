
context("est_map")

test_that("est_map for backcross autosome matches R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[1:19,]
    newmap <- est.map(hyper, err=0.002)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:19) {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("bc", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.002, 10000, 1e-6, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})

test_that("est_map for intercross autosome matches R/qtl", {

    library(qtl)
    data(listeria)
    listeria <- listeria[1:19,]
    newmap <- est.map(listeria, err=0.01)

    rf_qtl <- lapply(newmap, function(a) as.numeric(mf.h(diff(a))))

    for(i in 1:19) {
        g <- listeria$geno[[i]]$data
        g[is.na(g)] <- 0
        g <- g[rowSums(g!=0) > 1,] # omit individuals with <2 genotypes
        map <- listeria$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("f2", t(g), FALSE, rep(FALSE, nrow(g)),
                           matrix(ncol=nrow(g), nrow=0),
                           rf, 0.01, 10000, 1e-6, FALSE)

        expect_equivalent(rf_qtl[[i]], rf_qtl2)
        expect_equal(attr(newmap[[i]], "loglik"), attr(rf_qtl2, "loglik"))
   }

})


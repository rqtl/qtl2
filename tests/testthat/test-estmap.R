
context("est_map")

test_that("est_map for backcross autosome matches R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[4,]
    newmap <- est.map(hyper, err=0.002, maxit=100, tol=1e-6)
    rf_qtl <- as.numeric(mf.h(diff(newmap[[1]])))

    for(i in 1:1) {
        g <- hyper$geno[[i]]$data
        g[is.na(g)] <- 0
        map <- hyper$geno[[i]]$map
        rf <- mf.h(diff(map))

        rf_qtl2 <- est_map("bc", t(g), FALSE, rep(FALSE, nind(hyper)),
                           matrix(ncol=nind(hyper), nrow=0),
                           rf, 0.002, 100, 1e-6, TRUE)

        expect_equivalent(rf_qtl, rf_qtl2)
   }

})

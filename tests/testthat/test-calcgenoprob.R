
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


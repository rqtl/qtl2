context("scan1 with binary phenotype")

test_that("scan1 with binary phenotype gives same result as R/qtl", {

    library(qtl)
    data(hyper)
    hyper <- hyper[c(2,3,8),] # subset to 3 chromosomes
    hyper <- calc.genoprob(hyper, error.prob=0.002, step=1)
    hyper$pheno <- cbind(hyper$pheno, bp_binary=as.numeric(hyper$pheno[,1] > median(hyper$pheno[,1])))
    out1 <- scanone(hyper, pheno.col=3, model="binary", method="hk")

    # out1 map to list
    map <- split(out1[,2], out1[,1])
    mapn <- split(rownames(out1), out1[,1])
    for(i in seq_along(map)) names(map[[i]]) <- mapn[[i]]

    # convert for R/qtl2
    library(qtl2geno)
    hyper2 <- convert2cross2(hyper)
    pr <- calc_genoprob(hyper2, map, error_prob=0.002)
    out2 <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary")

    expect_equal(out1[,3] , as.numeric(out2))

    # add a couple of random covariates
    X <- matrix(rnorm(nind(hyper)*2, 20, 2), nrow=nind(hyper), ncol=2)
    rownames(X) <- ind_ids(hyper2)

    out1_ac <- scanone(hyper, pheno.col=3, model="binary", method="hk", addcovar=X)
    out2_ac <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary", addcovar=X)
    expect_equal(out1_ac[,3] , as.numeric(out2_ac))
    expect_equal(out1[,3] , as.numeric(out2))

})

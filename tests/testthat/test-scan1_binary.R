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
    hyper2 <- convert2cross2(hyper)
    pr <- calc_genoprob(hyper2, map, error_prob=0.002)
    out2 <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary")

    expect_equal(out1[,3] , as.numeric(out2))

    # add a couple of random covariates
    set.seed(45170702)
    X <- matrix(rnorm(nind(hyper)*2, 20, 2), nrow=nind(hyper), ncol=2)
    rownames(X) <- ind_ids(hyper2)

    out1_ac <- scanone(hyper, pheno.col=3, model="binary", method="hk", addcovar=X)
    out2_ac <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary", addcovar=X)
    expect_equal(out1_ac[,3] , as.numeric(out2_ac))

    # make one of them an interactive covariate
    out1_ic <- scanone(hyper, pheno.col=3, model="binary", method="hk",
                       addcovar=X, intcovar=X[,1])
    out2_ic <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary",
                     addcovar=X, intcovar=X[,1,drop=FALSE])
    expect_equal(out1_ic[,3] , as.numeric(out2_ic))

    # make both interactive covariates
    out1_ic2 <- scanone(hyper, pheno.col=3, model="binary", method="hk",
                        addcovar=X, intcovar=X)
    out2_ic2 <- scan1(pr, hyper2$pheno[,"bp_binary",drop=FALSE], model="binary",
                      addcovar=X, intcovar=X)
    expect_equal(out1_ic2[,3] , as.numeric(out2_ic2))

})


test_that("scan1 with binary phenotype in intercross gives same results as R/qtl", {

    skip_on_cran()

    library(qtl)
    data(listeria)
    listeria <- listeria[c(4,11,16),] # subset to 3 chromosomes
    listeria <- calc.genoprob(listeria, error.prob=0.002, step=1)
    listeria$pheno <- cbind(listeria$pheno, binary=as.numeric(listeria$pheno[,1] == 264))
    expect_warning(out1 <- scanone(listeria, pheno.col=3, model="binary", method="hk"))

    # out1 map to list
    map <- split(out1[,2], out1[,1])
    mapn <- split(rownames(out1), out1[,1])
    for(i in seq_along(map)) names(map[[i]]) <- mapn[[i]]

    # convert for R/qtl2
    listeria2 <- convert2cross2(listeria)
    pr <- calc_genoprob(listeria2, map, error_prob=0.002)
    out2 <- scan1(pr, listeria2$pheno[,"binary",drop=FALSE], model="binary")

    expect_equal(out1[,3] , as.numeric(out2))

    # add a couple of random covariates
    set.seed(45170702)
    X <- matrix(rnorm(nind(listeria)*2, 20, 2), nrow=nind(listeria), ncol=2)
    rownames(X) <- ind_ids(listeria2)

    expect_warning(out1_ac <- scanone(listeria, pheno.col=3, model="binary", method="hk", addcovar=X))
    out2_ac <- scan1(pr, listeria2$pheno[,"binary",drop=FALSE], model="binary", addcovar=X)
    expect_equal(out1_ac[,3] , as.numeric(out2_ac))

    # make one of them an interactive covariate
    expect_warning(out1_ic <- scanone(listeria, pheno.col=3, model="binary", method="hk",
                                      addcovar=X, intcovar=X[,1]))
    out2_ic <- scan1(pr, listeria2$pheno[,"binary",drop=FALSE], model="binary",
                     addcovar=X, intcovar=X[,1,drop=FALSE])
    expect_equal(out1_ic[,3] , as.numeric(out2_ic))

    # make both interactive covariates
    expect_warning(out1_ic2 <- scanone(listeria, pheno.col=3, model="binary", method="hk",
                                       addcovar=X, intcovar=X))
    out2_ic2 <- scan1(pr, listeria2$pheno[,"binary",drop=FALSE], model="binary",
                      addcovar=X, intcovar=X)
    expect_equal(out1_ic2[,3] , as.numeric(out2_ic2))

})


test_that("scan1 with multiple binary phenotype gives same results as R/qtl", {

    skip_on_cran()

    library(qtl)
    data(listeria)
    listeria <- listeria[c(1,3,12),] # subset to 3 chromosomes
    listeria <- calc.genoprob(listeria, error.prob=0.002, step=1)
    listeria$pheno <- cbind(listeria$pheno, binary1=as.numeric(listeria$pheno[,1] == 264),
                            binary2 = as.numeric(listeria$pheno[,1] > 116.5),
                            binary3 = as.numeric(listeria$pheno[,1] > 91.1))

    expect_warning(out1 <- scanone(listeria, pheno.col=3:5, model="binary", method="hk"))

    # out1 map to list
    map <- split(out1[,2], out1[,1])
    mapn <- split(rownames(out1), out1[,1])
    for(i in seq_along(map)) names(map[[i]]) <- mapn[[i]]

    # convert for R/qtl2
    listeria2 <- convert2cross2(listeria)
    pr <- calc_genoprob(listeria2, map, error_prob=0.002)
    out2 <- scan1(pr, listeria2$pheno[,2:4], model="binary")

    for(i in 1:3) expect_equal(out1[,i+2], as.numeric(out2[,i]))

    # add a couple of random covariates
    set.seed(45170702)
    X <- matrix(rnorm(nind(listeria)*2, 20, 2), nrow=nind(listeria), ncol=2)
    rownames(X) <- ind_ids(listeria2)

    expect_warning(out1_ac <- scanone(listeria, pheno.col=3:5, model="binary", method="hk", addcovar=X))
    out2_ac <- scan1(pr, listeria2$pheno[,2:4], model="binary", addcovar=X)
    for(i in 1:3) expect_equal(out1_ac[,i+2], as.numeric(out2_ac[,i]))

    # make one of them an interactive covariate
    expect_warning(out1_ic <- scanone(listeria, pheno.col=3:5, model="binary", method="hk",
                                      addcovar=X, intcovar=X[,1]))
    out2_ic <- scan1(pr, listeria2$pheno[,2:4], model="binary",
                     addcovar=X, intcovar=X[,1,drop=FALSE])
    for(i in 1:3) expect_equal(out1_ic[,i+2], as.numeric(out2_ic[,i]))

    # make both interactive covariates
    expect_warning(out1_ic2 <- scanone(listeria, pheno.col=3:5, model="binary", method="hk",
                                       addcovar=X, intcovar=X))
    out2_ic2 <- scan1(pr, listeria2$pheno[,2:4], model="binary",
                      addcovar=X, intcovar=X)
    for(i in 1:3) expect_equal(out1_ic2[,i+2], as.numeric(out2_ic2[,i]))

})


test_that("scan1 with weights gives the same answer as glm()", {

    skip_on_cran()

    library(qtl)
    data(listeria)
    set.seed(20180717)
    listeria$pheno$sex <- sample(0:1, nind(listeria), replace=TRUE)
    listeria <- listeria[c(1,3,"X"), ] # subset to 3 chromosomes
    listeria <- convert2cross2(listeria)
    Xcovar <- get_x_covar(listeria)

    phe <- cbind(binary1=as.numeric(listeria$pheno[,1] == 264),
                 binary2 = as.numeric(listeria$pheno[,1] > 116.5),
                 binary3 = as.numeric(listeria$pheno[,1] > 91.1))
    rownames(phe) <- rownames(listeria$pheno)

    map <- insert_pseudomarkers(listeria$gmap, step=2.5)
    pr <- calc_genoprob(listeria, map)

    out <- scan1(pr, phe, model="binary", Xcovar=Xcovar)

    nmar <- sapply(map, length) # no. markers/pseudomarkers per chromosome
    csmar <- cumsum(c(0, nmar))
    pos <- c(37, 7, 9)

    # null loglik (X chr needs sex as covariate)
    dev_glm0 <- apply(phe, 2, function(y)
        glm(y ~ 1, family=binomial(link=logit))$deviance)
    dev_glm0X <- apply(phe, 2, function(y)
        glm(y ~ Xcovar, family=binomial(link=logit))$deviance)
    dev_glm0 <- list(dev_glm0, dev_glm0, dev_glm0X)

    for(i in 1:3) {
        dev_glm <- apply(phe, 2, function(y)
            glm(y ~ -1 + pr[[i]][,,pos[i]], family=binomial(link=logit))$deviance)

        lod_glm <- (dev_glm0[[i]] - dev_glm)/(2*log(10))

        expect_equal(out[csmar[i] + pos[i],], lod_glm)
    }

    # repeat with weights
    weights <- setNames(sample(1:4, n_ind(listeria), replace=TRUE), rownames(phe))

    out <- scan1(pr, phe, model="binary", Xcovar=Xcovar, weights=weights)

    # null loglik (X chr needs sex as covariate)
    dev_glm0 <- apply(phe, 2, function(y)
        glm(y ~ 1, family=binomial(link=logit), weights=weights)$deviance)
    dev_glm0X <- apply(phe, 2, function(y)
        glm(y ~ Xcovar, family=binomial(link=logit), weights=weights)$deviance)
    dev_glm0 <- list(dev_glm0, dev_glm0, dev_glm0X)

    for(i in 1:3) {
        dev_glm <- apply(phe, 2, function(y)
            glm(y ~ -1 + pr[[i]][,,pos[i]], family=binomial(link=logit), weights=weights)$deviance)

        lod_glm <- (dev_glm0[[i]] - dev_glm)/(2*log(10))

        expect_equal(out[csmar[i] + pos[i],], lod_glm)
    }

})

context("reduced rank covariates")

test_that("scan1 etc work with reduced-rank covariates", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(1,12,"X")]

    Xcovar <- get_x_covar(iron)
    set.seed(3706216)
    n_batch <- 4
    batch <- sample(1:n_batch, n_ind(iron), replace=TRUE)
    X <- matrix(0, nrow=n_ind(iron), ncol=n_batch-1)
    dimnames(X) <- list(rownames(Xcovar), paste0("batch", 2:n_batch))
    for(i in 2:n_batch) X[batch==i,i-1] <- 1

    pr <- calc_genoprob(iron, error_prob=0.01)

    phe <- iron$pheno
    phe[batch==n_batch,1] <- NA
    phe[batch==1,2] <- NA

    k <- calc_kinship(pr)

    # scan1 with no kinship
    out <- scan1(pr, phe, Xcovar=Xcovar, addcovar=X)
    expected <- cbind(scan1(pr, phe[!is.na(phe[,1]),1,drop=FALSE], Xcovar=Xcovar, addcovar=X[,-3]),
                      scan1(pr, phe[!is.na(phe[,2]),2,drop=FALSE], Xcovar=Xcovar, addcovar=X[,-3]))
    expect_equal(out, expected)

    # scan1 with kinship
    out_k <- scan1(pr, phe, k, Xcovar=Xcovar, addcovar=X)
    expected_k <- cbind(scan1(pr, phe[!is.na(phe[,1]),1,drop=FALSE], k, Xcovar=Xcovar, addcovar=X[,-3]),
                        scan1(pr, phe[!is.na(phe[,2]),2,drop=FALSE], k, Xcovar=Xcovar, addcovar=X[,-3]))
    expect_equal(out_k, expected_k, tol=5e-6)

    # scan1coef with no kinship
    for(phecol in 1:2) {
        for(chr in names(pr)) {
            co <- scan1coef(pr[,chr], phe[,phecol,drop=FALSE], addcovar=X)
            Xcol2keep <- colnames(X)[colnames(X) %in% colnames(co)]
            expected <- scan1coef(pr[,chr], phe[!is.na(phe[,phecol]),phecol,drop=FALSE], 
                                  addcovar=X[,Xcol2keep])
            expect_equal(co, expected)
        }
    }

    # scan1coef with kinship
    for(phecol in 1:2) {
        for(chr in names(pr)) {
            co <- scan1coef(pr[,chr], phe[,phecol,drop=FALSE], k, addcovar=X)
            Xcol2keep <- colnames(X)[colnames(X) %in% colnames(co)]
            expected <- scan1coef(pr[,chr], phe[!is.na(phe[,phecol]),phecol,drop=FALSE], 
                                  k, addcovar=X[,Xcol2keep])
            expect_equal(co, expected)
        }
    }


    # scan1blup with no kinship
    for(phecol in 1:2) {
        for(chr in names(pr)) {
            co <- scan1blup(pr[,chr], phe[,phecol,drop=FALSE], addcovar=X)
            Xcol2keep <- colnames(X)[colnames(X) %in% colnames(co)]
            expected <- scan1blup(pr[,chr], phe[!is.na(phe[,phecol]),phecol,drop=FALSE], 
                                  addcovar=X[,Xcol2keep])
            expect_equal(co, expected)
        }
    }

    # scan1blup with kinship
    for(phecol in 1:2) {
        for(chr in names(pr)) {
            co <- scan1blup(pr[,chr], phe[,phecol,drop=FALSE], k, addcovar=X)
            Xcol2keep <- colnames(X)[colnames(X) %in% colnames(co)]
            expected <- scan1blup(pr[,chr], phe[!is.na(phe[,phecol]),phecol,drop=FALSE], 
                                  k, addcovar=X[,Xcol2keep])
            expect_equal(co, expected)
        }
    }

})

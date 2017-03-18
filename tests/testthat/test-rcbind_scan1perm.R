context("rbind/cbind scan1perm")

# load qtl2geno package for data and genoprob calculation
library(qtl2geno)

# read data
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
iron <- iron[,c(18,19,"X")]

# insert pseudomarkers into map
map <- insert_pseudomarkers(iron$gmap, step=1)

# calculate genotype probabilities
probs <- calc_genoprob(iron, map, error_prob=0.002)

# grab phenotypes and covariates; ensure that covariates have names attribute
pheno <- iron$pheno
covar <- match(iron$covar$sex, c("f", "m")) # make numeric
names(covar) <- rownames(iron$covar)
Xcovar <- get_x_covar(iron)

# strata for permutations
perm_strata <- mat2strata(Xcovar)

test_that("rbind.scan1perm works", {

    # not X-chr-specific
    operm1 <- scan1perm(probs, pheno, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    operm2 <- scan1perm(probs, pheno, Xcovar=Xcovar, n_perm=4, perm_strata=perm_strata)
    operm <- rbind(operm1, operm2)
    expected <- rbind(unclass(operm1), unclass(operm2))
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with three
    operm <- rbind(operm1, operm2, operm1)
    expected <- rbind(unclass(operm1), unclass(operm2), unclass(operm1))
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with one
    expect_equal(rbind(operm1), operm1)

    # X-chr-specific
    operm1 <- scan1perm(probs, pheno, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    operm2 <- scan1perm(probs, pheno, Xcovar=Xcovar, n_perm=4, perm_strata=perm_strata,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    operm <- rbind(operm1, operm2)
    expected <- list(A=rbind(operm1$A, operm2$A),
                     X=rbind(operm1$X, operm2$X))
    attr(expected, "chr_lengths") <- attr(operm1, "chr_lengths")
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with three
    operm <- rbind(operm1, operm2, operm1)
    expected <- list(A=rbind(operm1$A, operm2$A, operm1$A),
                     X=rbind(operm1$X, operm2$X, operm1$X))
    attr(expected, "chr_lengths") <- attr(operm1, "chr_lengths")
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with one
    expect_equal(rbind(operm1), operm1)

})


test_that("cbind.scan1perm works", {

    # not X-chr-specific, same number of permutations
    operm1 <- scan1perm(probs, pheno[,1,drop=FALSE], Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    operm2 <- scan1perm(probs, pheno[,2,drop=FALSE], Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    operm <- cbind(operm1, operm2)
    expected <- cbind(unclass(operm1), unclass(operm2))
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with three
    operm <- cbind(operm1, operm2, operm1)
    expected <- cbind(unclass(operm1), unclass(operm2), unclass(operm1))
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with one
    expect_equal(cbind(operm1), operm1)



    # different numbers of permutations
    operm2b <- operm2[1:2,,drop=FALSE]
    operm <- cbind(operm1, operm2b)
    expected <- cbind(operm1, operm2)
    expected[3,2] <- NA
    expect_equal(operm, expected)



    # X-chr-specific
    operm1 <- scan1perm(probs, pheno[,1,drop=FALSE], Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    operm2 <- scan1perm(probs, pheno[,2,drop=FALSE], Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    operm <- cbind(operm1, operm2)
    expected <- list(A=cbind(operm1$A, operm2$A),
                     X=cbind(operm1$X, operm2$X))
    attr(expected, "chr_lengths") <- attr(operm1, "chr_lengths")
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with three
    operm <- cbind(operm1, operm2, operm1)
    expected <- list(A=cbind(operm1$A, operm2$A, operm1$A),
                     X=cbind(operm1$X, operm2$X, operm1$X))
    attr(expected, "chr_lengths") <- attr(operm1, "chr_lengths")
    class(expected) <- class(operm1)
    expect_equal(operm, expected)

    # test with one
    expect_equal(cbind(operm1), operm1)



    # different numbers of permutations
    operm2b <- list(A=operm2$A[1:2,,drop=FALSE],
                    X=operm2$X[1:5,,drop=FALSE])
    attr(operm2b, "chr_lengths") <- attr(operm2, "chr_lengths")
    operm <- cbind(operm1, operm2b)
    expected <- cbind(operm1, operm2)
    expected$A[3,2] <- NA
    expected$X[6:7,2] <- NA
    expect_equal(operm, expected)

})

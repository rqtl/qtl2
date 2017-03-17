context("summary scan1perm")

test_that("summary_scan1perm works", {

    seed <- 9896433
    RNGkind("Mersenne-Twister")

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

    # not X-chr-specific
    set.seed(seed)
    operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar, n_perm=3)

    result <- summary(operm)
    expected <- rbind("0.05" = c(liver=1.5480532515070430932,spleen=1.1945482736887691466))
    attr(expected, "n_perm") <- rbind(c(liver=3, spleen=3))
    class(expected) <- c("summary.scan1perm", "matrix")
    expect_equal(result, expected)

    # two thresholds
    result <- summary(operm, c(0.05, 0.2))
    expected <- rbind(expected, "0.2"=c(liver=1.5039885621071871213, spleen=1.0497165924333802245))
    attr(expected, "n_perm") <- rbind(c(liver=3, spleen=3))
    class(expected) <- c("summary.scan1perm", "matrix")
    expect_equal(result, expected)

    # X-chr-specific
    set.seed(seed)
    operm <- scan1perm(probs, pheno, addcovar=covar, Xcovar=Xcovar,
                       n_perm=3, perm_Xsp=TRUE,
                       chr_lengths=chr_lengths(iron$gmap))
    result <- summary(operm)
    expected <- list(A=rbind("0.05" = c(liver=1.5244058138572063044, spleen=0.63833353205793441631)),
                     X=rbind("0.05" = c(liver=1.7400163974820828106, spleen= 2.8361536028418536937)))
    attr(expected, "n_perm") <- rbind(A=c(liver=3, spleen=3), X=c(liver=7, spleen=7))
    class(expected) <- c("summary.scan1perm", "list")
    expect_equal(result, expected)

    # two thresholds
    result <- summary(operm, c(0.05, 0.2))
    expected$A <- rbind(expected$A, "0.2"=c(liver=1.4053300973638913618, spleen=0.61906591053182846718))
    expected$X <- rbind(expected$X, "0.2"=c(liver=1.6152559270167379246, spleen=2.2975884236147399164))
    expect_equal(result, expected)

})

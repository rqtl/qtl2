context("fit1")

test_that("fit1 by H-K works in intercross", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(18:19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)
    Xcovar <- get_x_covar(iron)

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov) })

    # fit1, no missing data
    npos <- sapply(probs, function(a) dim(a)[3])
    pmar <- c(3, 4, 12)
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov) })

    pos <- cumsum(c(0, npos[-3])) + pmar
    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:3) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:3)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])

    # repeat the whole thing with a couple of missing phenotypes
    pheno[c(187, 244)] <- NA

    # calculate LOD scores
    out <- scan1(probs, pheno, addcovar=covar, Xcovar=Xcovar)

    # estimate coefficients; no covariates for X chromosome
    coef <- lapply(seq_len(length(probs)), function(i) {
        if(i==3) cov <- NULL
        else cov <- covar
        scan1coef(subset(probs, chr=names(probs)[i]), pheno, addcovar=cov) })

    # fit1, no missing data
    out_fit1 <- lapply(seq(along=pmar),
                       function(i) {
        if(i==3) { nullcov <- Xcovar; cov <- NULL } # need Xcovar under null on X chr but no other covariates
        else { nullcov <- NULL; cov <- covar }      # sex as covariate; no additional covariates under null
        fit1(probs[[i]][,,pmar[i]], pheno, addcovar=cov, nullcovar=nullcov) })

    # check LOD vs scan1, plus ind'l contributions to LOD
    for(i in 1:3) {
        expect_equal(out_fit1[[i]]$lod, out[pos[i],1])
        expect_equal(sum(out_fit1[[i]]$ind_lod), out_fit1[[i]]$lod)
    }

    # check coefficients
    for(i in 1:3)
        expect_equal(out_fit1[[i]]$coef, coef[[i]][pmar[i],])
})

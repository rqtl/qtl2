context("scan1gen")

test_that("scan1gen works", {

    # function to fit glm
    ll_glm <-
        function(pr, pheno, addcovar=NULL, ...)
        {
            if(ncol(pheno)>1)
                return(sapply(1:ncol(pheno), function(i) ll_glm(pr, pheno[,i,drop=FALSE], addcovar=addcovar, ...)))

            formula <- ifelse(is.null(pr), "pheno ~ 1", "pheno ~ pr")
            if(!is.null(addcovar)) formula <- paste(formula, "+ addcovar")

            glm_out <- glm(as.formula(formula), family=binomial(link=logit))
            -glm_out$deviance/(2*log(10)) # log10 likelihood
        }

    # read data
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c("1", "16", "X")] # subset to chr 1, 16, and X

    # insert pseudomarkers into map
    map <- insert_pseudomarkers(iron$gmap, step=1)

    # calculate genotype probabilities
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    # covariates for X chr under null
    Xcovar <- get_x_covar(iron)

    # create binary trait
    bin_pheno <- setNames(as.numeric(iron$pheno[,1] > median(iron$pheno[,1])),
                          rownames(iron$pheno))

    expect_equal(scan1gen(probs, bin_pheno, Xcovar=Xcovar, func=ll_glm),
                 scan1(probs, bin_pheno, Xcovar=Xcovar, model="binary"))


    # try a second binary trait
    bin_pheno <- cbind(bin_pheno,
                       as.numeric(iron$pheno[,1] > quantile(iron$pheno[,1], 0.75)))

    # vectorize the function
    expect_equal(scan1gen(probs, bin_pheno, Xcovar=Xcovar, func=ll_glm),
                 scan1(probs, bin_pheno, Xcovar=Xcovar, model="binary"))

    # don't vectorize the function
    expect_equal(scan1gen(probs, bin_pheno, Xcovar=Xcovar, func=ll_glm, vectorize_func=FALSE),
                 scan1(probs, bin_pheno, Xcovar=Xcovar, model="binary"))

})

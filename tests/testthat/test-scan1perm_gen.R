context("scan1 permutations with scan_func")

test_that("scan1 permutations work with scan_func", {

    skip_if(isnt_karl(), "this test only run locally")

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    Xcovar <- get_x_covar(iron)
    bin_pheno <- setNames(as.numeric(iron$pheno[,1] > median(iron$pheno[,1])),
                          rownames(iron$pheno))

    scanf <- function(genoprobs, pheno, kinship, addcovar, Xcovar, intcovar, weights)
        scan1(genoprobs, pheno=pheno, kinship=kinship, addcovar=addcovar,
              Xcovar=Xcovar, intcovar=intcovar, weights=weights, model="binary")

    set.seed(20260615)
    operm1 <- scan1perm(probs, bin_pheno, Xcovar=Xcovar, model="binary", n_perm=10)
    set.seed(20260615)
    operm2 <- scan1perm(probs, bin_pheno, Xcovar=Xcovar, scan_func=scanf, n_perm=10)

    expect_equal(operm1, operm2)

    ### glm with link=probit

    ll_glm <-
    function(pr, pheno, addcovar=NULL, ...)
    {
        formula <- ifelse(is.null(pr), "pheno ~ 1", "pheno ~ pr")
        if(!is.null(addcovar)) formula <- paste(formula, "+ addcovar")

        glm_out <- glm(as.formula(formula), family=binomial(link=probit))
        -glm_out$deviance/(2*log(10)) # log10 likelihood
    }

    set.seed(20260615)
    operm3 <- scan1perm(probs, bin_pheno, Xcovar=Xcovar, scan_func=scan1gen,
                        func=ll_glm, n_perm=10, cores=10)

    expected <- structure(c(2.17385652559186, 2.31427078248578, 3.43941883075001,
                            1.91898059153991, 2.01118425059423, 2.32202927408811, 1.72215414380619,
                            1.31001918257803, 1.37739899949347, 2.68954157272938),
                          dim = c(10L, 1L), dimnames = list(NULL, "pheno1"),
                          class = c("scan1perm", "matrix"))
    expect_equal(operm3, expected)

    set.seed(20260615)
    operm4 <- scan1perm(probs[,1:19], bin_pheno, scan_func=scan1gen,
                        func=ll_glm, n_perm=10, cores=10)

    expected <- structure(c(2.68609128205085, 1.77002204352716, 1.93368949267834,
                            2.31094501376531, 2.79567905567447, 1.79725590191482, 2.79244070919604,
                            1.3401617230496, 3.955798656798, 2.23867940226528),
                          dim = c(10L, 1L), dimnames = list(NULL, "pheno1"),
                          class = c("scan1perm", "matrix"))
    expect_equal(operm4, expected)


    set.seed(20260615)
    operm5 <- scan1perm(probs, bin_pheno, Xcovar=Xcovar, scan_func=scan1gen,
                        func=ll_glm, n_perm=10, cores=10, perm_Xsp=TRUE,
                        chr_lengths=chr_lengths(iron$gmap))

    expected_summary <- structure(list(A = structure(c(2.0702243751226, 2.20882194904352, 2.56447872617982, 2.60840127327936),
                                                     dim = c(4L, 1L), dimnames = list(c("0.5", "0.25", "0.1", "0.05"), "pheno1")),
                                       X = structure(c(2.79707347743163, 3.33791604605526, 3.61256992679585, 4.07772316813587),
                                                     dim = c(4L, 1L), dimnames = list(c("0.5", "0.25", "0.1", "0.05"), "pheno1"))),
                                  class = c("summary.scan1perm", "list"),
                                  n_perm = structure(c(10, 283), dim = 2:1, dimnames = list(c("A", "X"), "pheno1")))
    expect_equal(summary(operm5, c(0.5, 0.25, 0.1, 0.05)), expected_summary)


    # permutations with scan1snps
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    DOex <- DOex[,c("2", "3")] # subset to chr 2 and 3
    probs <- calc_genoprob(DOex, error_prob=0.002)

    snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
    queryf <- create_variant_query_func(snpdb_file)

    set.seed(20260616)
    operm6 <- scan1perm(genoprobs=probs, map=DOex$pmap, pheno=DOex$pheno,
                       scan_func=scan1snps, query_func=queryf, n_perm=10)
    expected <- structure(c(1.87603485437823, 2.78903657330818, 1.28168354734608,
                            1.96766231770538, 2.00872833761448, 2.64571234938997, 2.84375635883363,
                            1.86947672509228, 1.49270250736036, 1.7635254846055),
                          dim = c(10L, 1L), dimnames = list(NULL, "OF_immobile_pct"),
                          class = c("scan1perm", "matrix"))
    expect_equal(operm6, expected)

})

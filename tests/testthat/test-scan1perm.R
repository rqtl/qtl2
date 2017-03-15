context("scan1 permutations")

test_that("scan1 permutations work (regression test)", {

    seed <- 3025685

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata","iron.zip", package="qtl2geno"))
    iron <- iron[,c(18,19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, map, err=0.002)

    pheno1 <- iron$pheno
    pheno2 <- iron$pheno; pheno2[5,1] <- NA

    Xcovar <- get_x_covar(iron)
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)
    perm_strat <- apply(Xcovar, 1, paste, collapse=":")

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, n_perm=3)
    expected <- cbind(liver= c(1.3756640664429049536, 2.7572772758036521168, 0.52452198541287131661),
                      spleen=c(1.5153608426190210423, 1.4532027591584668613, 2.86857866542299611012))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strat=perm_strat)
    expected <- cbind(liver= c(1.56882977009222202014, 1.0308336167670439920, 1.6771296752179409850),
                      spleen=c(0.88478150549593514995, 1.1297760401050815915, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strat=perm_strat)
    expected <- cbind(liver= c(2.0685972274435684426, 1.5683121279062959275, 1.6771296752180671064),
                      spleen=c(1.0328759662065980507, 1.5091023963920129347, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strat=perm_strat)
    expected <- list(A=cbind(liver= c(1.56882977009222202014, 0.17654731701161452406, 0.67967768984640919427),
                             spleen=c(0.88478150549593514995, 1.12977604010508159149, 0.95835278347981578406)),
                     X=cbind(liver= c(0.65635525818365536566, 0.59711143944156575003, 0.49516910044845729999, 0.70305870212073351411,
                                      0.44456114810071056809, 1.17745924272541380160, 0.69147359544272646303),
                             spleen=c(0.71360878105542902006, 0.54087014463646632123, 1.56708158211767667467, 1.21950357033675516050,
                                      1.49396169005535739416, 2.8719171298632524270, 1.08820611768204855707)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)

    # one missing phenotype, no covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, n_perm=3)
    expected <- cbind(liver= c(1.4246052689684201020, 1.8035976886303082267, 1.295827488548965345),
                      spleen=c(2.1450871714236026122, 1.9160146558007777884, 1.887126247743548646))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strat=perm_strat)
    expected <- cbind(liver= c(1.56455340995498248802, 1.0121930261312361843, 1.6731999914352853054),
                      spleen=c(0.88478150549593514995, 1.1297760401050815915, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strat=perm_strat)
    expected <- cbind(liver= c(2.0627670956684416304, 1.6261084401604977145, 1.6731999914352853054),
                      spleen=c(1.0328759662065980507, 1.5091023963920129347, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    expected <- list(A=cbind(liver= c(1.18430970191987361420, 0.85274796222984994287,  0.76981903233409942899),
                             spleen=c(1.96074450827722834840, 2.47749749792303752827,  1.20235911795879424346)),
                     X=cbind(liver= c(1.29867858019373416670, 0.25258741397366168968,  2.22594902965746355150, 0.69668535735455083824,
                                      1.84525108476481491948, 1.97947969062733486467, -0.07964898968422540193),
                             spleen=c(1.39377789098622173470, 1.24497337934968221873,  0.94154476905226225369, 0.50568976179249958136,
                                      0.70117350389718069437, 0.46582839190887348479,  0.98324516512323434370)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)

})

context("scan1 permutations")

test_that("scan1 permutations work (regression test)", {

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

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
    perm_strata <- mat2strata(Xcovar)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, n_perm=3)
    expected <- cbind(liver= c(1.3756640664429049536, 2.7572772758036521168, 0.52452198541287131661),
                      spleen=c(1.5153608426190210423, 1.4532027591584668613, 2.86857866542299611012))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, n_perm=3, perm_strata=perm_strata)
    expected <- cbind(liver= c(13.541670362333418254, 14.7024286141676103767, 14.1453256570137018144),
                      spleen=c( 5.930376703877245248,  7.4039049500925262493,  5.5851374704542742222))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.56882977009222202014, 1.0308336167670439920, 1.6771296752179409850),
                      spleen=c(0.88478150549593514995, 1.1297760401050815915, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(2.0685972274435684426, 1.5683121279062959275, 1.6771296752180671064),
                      spleen=c(1.0328759662065980507, 1.5091023963920129347, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
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
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strata=perm_strata)
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
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.56455340995498248802, 1.0121930261312361843, 1.6731999914352853054),
                      spleen=c(0.88478150549593514995, 1.1297760401050815915, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(2.0627670956684416304, 1.6261084401604977145, 1.6731999914352853054),
                      spleen=c(1.0328759662065980507, 1.5091023963920129347, 2.8456095237703333822))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    perm_strata <- rep(1, nrow(pheno2)) # avoid the stratified permutations
    names(perm_strata) <- rownames(pheno2)
    operm <- scan1perm(pr, pheno2, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                       perm_Xsp=TRUE, chr_lengths=chr_lengths(map),
                       perm_strata=perm_strata)
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


test_that("scan1 permutations work with single kinship matrix (regression test)", {

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata","iron.zip", package="qtl2geno"))
    iron <- iron[,c(18,19,"X")]
    map <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, map, err=0.002)
    kinship <- calc_kinship(pr)

    pheno1 <- iron$pheno
    pheno2 <- iron$pheno; pheno2[5,1] <- NA

    Xcovar <- get_x_covar(iron)
    sex <- as.numeric(iron$covar$sex=="m")
    names(sex) <- rownames(iron$covar)
    perm_strata <- mat2strata(Xcovar)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, n_perm=3)
    expected <- cbind(liver= c(1.5850327051161108294, 1.6225465841443282855, 1.1047230815073763033),
                      spleen=c(2.1929283392327243440, 1.9397727636298436327, 1.8698391778175107447))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, n_perm=3, perm_strata=perm_strata)
    expected <- cbind(liver= c(14.6442198530459659622, 14.6817005229127524046, 15.6345815802552241536),
                      spleen=c( 5.5167589745319025596,  5.6328648432420660441,  8.0494249754664739527))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.29766020287924632726, 1.0197730259964858934, 2.0112145151147022837),
                      spleen=c(0.90427196894971517693, 1.0214860518121018362, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.8645265744331560587, 1.6155547860418919548, 2.0112145151147022837),
                      spleen=c(1.0356350517366044173, 1.4490075143563088123, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    expected <- list(A=cbind(liver= c(1.29766020287924632726, 0.14126238719321179693, 0.77580567004437706036),
                             spleen=c(0.90427196894971517693, 1.02148605181210183623, 0.93168604007222088903)),
                     X=cbind(liver= c(0.73307723185483875117, 0.99760458934328322123, 0.83443976911742778757, 0.89370279069466052047,
                                      0.67087690411915201771, 1.50033439242864674590, 0.88505044142934741203),
                             spleen=c(0.83259582111266161597, 0.61786570853812217141, 1.62939070897649207481, 1.30371105730912195675,
                                      1.58674017274476608641, 2.91485391258331727470, 1.16090441720256909441)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, no covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, n_perm=3)
    expected <- cbind(liver= c(1.5863339748532325757, 1.6385708572073809375, 1.1068373425454520742),
                      spleen=c(2.1929283392327243440, 1.9397727636298436327, 1.8698391778175107447))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.29518706542233608126, 0.99913640281598758985, 2.0054677281986590387),
                      spleen=c(0.90427196894971517693, 1.02148605181210183623, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.8604662282096269266, 1.6877912266523966700, 2.0054677281986590387),
                      spleen=c(1.0356350517366044173, 1.4490075143563088123, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    perm_strata <- rep(1, nrow(pheno2)) # avoid the stratified permutations
    names(perm_strata) <- rownames(pheno2)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                       perm_Xsp=TRUE, chr_lengths=chr_lengths(map),
                       perm_strata=perm_strata)
    expected <- list(A=cbind(liver= c(1.37030487476505635770, 0.58183562848782044430, 0.79833486499453443219),
                             spleen=c(2.01478095022342573730, 2.52716001067115758620, 1.26929860408445516207)),
                     X=cbind(liver= c(1.87350146601092815150, 0.91524542868117597649, 2.3045977403837074604, 1.3042235679000910853,
                                      2.12568139604854344782, 2.37620592723849499706, 0.32068152868385429999),
                             spleen=c(1.68345573592034480900, 1.74381613053044959294, 1.23338578366108553300, 0.9426755474293885273,
                                      0.93450566468420337429, 0.89848750594932114133, 1.47549973926410826763)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)

})

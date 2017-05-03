context("scan1 permutations")

library(qtl2geno)
iron <- read_cross2(system.file("extdata","iron.zip", package="qtl2geno"))
iron <- iron[,c(18,19,"X")]
map <- insert_pseudomarkers(iron$gmap, step=1)
pr <- calc_genoprob(iron, map, err=0.002)
kinship <- calc_kinship(pr)
kinship_loco <- calc_kinship(pr, "loco")

pheno1 <- iron$pheno
pheno2 <- iron$pheno; pheno2[5,1] <- NA

Xcovar <- get_x_covar(iron)
sex <- as.numeric(iron$covar$sex=="m")
names(sex) <- rownames(iron$covar)
perm_strata <- mat2strata(Xcovar)

test_that("scan1 permutations work (regression test)", {

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

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

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, n_perm=3)
    expected <- cbind(liver= c(1.5850327051161108294, 1.6225465841443282855, 1.1047230815073763033),
                      spleen=c(2.1929283392327243440, 1.9397727636298436327, 1.8698391778175107447))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, n_perm=3, perm_strata=perm_strata)
    expected <- cbind(liver= c(14.6442198530459659622, 14.6817005229127524046, 15.6345815802552241536),
                      spleen=c( 5.5167589745319025596,  5.6328648432420660441,  8.0494249754664739527))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.29766020287924632726, 1.0197730259964858934, 2.0112145151147022837),
                      spleen=c(0.90427196894971517693, 1.0214860518121018362, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.8645265744331560587, 1.6155547860418919548, 2.0112145151147022837),
                      spleen=c(1.0356350517366044173, 1.4490075143563088123, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    expected <- list(A=cbind(liver= c(1.29766020287924632726, 0.14126238719321179693, 0.77580567004437706036),
                             spleen=c(0.90427196894971517693, 1.02148605181210183623, 0.93168604007222088903)),
                     X=cbind(liver= c(0.50034270818987047758, 0.76231900452964485027, 0.6011447003868983785, 0.65860596461113896094,
                                      0.43605696449176051255, 1.26565994163476069900, 0.65084718581631106904),
                             spleen=c(0.73143027950382133451, 0.53598997526824987414, 1.5569757625108506804, 1.21787150923781783973,
                                      1.49387241817333116245, 2.83703932517150425600, 1.07984190845146876825)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected, tolerance=3e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=3e-7)

    # one missing phenotype, no covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, n_perm=3)
    expected <- cbind(liver= c(1.5863339748532325757, 1.6385708572073809375, 1.1068373425454520742),
                      spleen=c(2.1929283392327243440, 1.9397727636298436327, 1.8698391778175107447))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)

    # one missing phenotype, sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.29518706542233608126, 0.99913640281598758985, 2.0054677281986590387),
                      spleen=c(0.90427196894971517693, 1.02148605181210183623, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=2e-7)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.8604662282096269266, 1.6877912266523966700, 2.0054677281986590387),
                      spleen=c(1.0356350517366044173, 1.4490075143563088123, 2.8832996501350738328))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=1e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=1e-7)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    perm_strata <- rep(1, nrow(pheno2)) # avoid the stratified permutations
    names(perm_strata) <- rownames(pheno2)
    operm <- scan1perm(pr, pheno2, kinship, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                       perm_Xsp=TRUE, chr_lengths=chr_lengths(map),
                       perm_strata=perm_strata)
    expected <- list(A=cbind(liver= c(1.37030487476505635770, 0.58183562848782044430, 0.79833486499453443219),
                             spleen=c(2.01478095022342573730, 2.52716001067115758620, 1.26929860408445516207)),
                     X=cbind(liver= c(1.63673834214325419900, 0.67953625118790006443, 2.073666671975550102000, 1.06998409582742026736,
                                      1.89434763183230137074, 2.14236446791721579785, 0.086724900384358594163),
                             spleen=c(1.56170484328174996590, 1.56288576386423483378, 1.116409713761464583800, 0.78539826876917206988,
                                      0.82004164028519510588, 0.74049399535694393482, 1.310035745047248845196)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected, tolerance=1e-7)

})


test_that("scan1 permutations work with LOCO kinship matrix (regression test)", {

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, n_perm=3)
    expected <- cbind(liver= c(1.4233905282897871825, 1.6883979222096348050, 1.1105319041474899233),
                      spleen=c(2.1849905121197314983, 1.9112972630703948251, 1.8871262477436099303))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=1e-7)

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, n_perm=3, perm_strata=perm_strata)
    expected <- cbind(liver= c(15.0498110236199185152, 15.1559669590417573914, 16.1619790013112449856),
                      spleen=c( 5.6528544325570964091,  5.7416707117225200818,  8.1270354686771071329))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.56882977009233215426, 1.0513515452269666106, 2.0908041328949908966),
                      spleen=c(0.91750143650706161846, 1.0096069617524108253, 2.8456095237817189414))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(2.0685972274436985607, 1.5532309407813880142, 2.0908041328949908966),
                      spleen=c(1.0453315747253293377, 1.4393283407079127123, 2.8456095237817189414))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=2e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=2e-7)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    expected <- list(A=cbind(liver= c(1.56882977009233215426, 0.15408564377888087082, 0.67967768984653986752),
                             spleen=c(0.91750143650706161846, 1.00960696175241082528, 0.91230958640896986367)),
                     X=cbind(liver= c(0.47525078870451631374, 0.81648974630841897326, 0.55021896067082298742, 0.68285545241312051168,
                                      0.46054636746004662395, 1.28223568383226194100, 0.66115351693360990826),
                             spleen=c(0.71360878106659053621, 0.54087014464766181021, 1.56708158212895898309, 1.21950357034804057754,
                                      1.49396169006666879042, 2.87191712987437774980, 1.08820611769328978724)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected, tolerance=3e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=3e-7)

    # one missing phenotype, no covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship_loco, n_perm=3)
    expected <- cbind(liver= c(1.4246052689684871595, 1.7056415036570125032, 1.1168979618114607266),
                      spleen=c(2.1849905121197314983, 1.9112972630703948251, 1.8871262477436099303))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=1e-7)

    # one missing phenotype, sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.56455340995486880118, 1.0297698518383779920, 2.0850243087458522062),
                      spleen=c(0.91750143650706161846, 1.0096069617524108253, 2.8456095237817189414))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=1e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=1e-7)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship_loco, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(2.0627670956682830905, 1.6203949273894016070, 2.0850243087458522062),
                      spleen=c(1.0453315747253293377, 1.4393283407079127123, 2.8456095237817189414))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected, tolerance=1e-7)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, kinship_loco, addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected, tolerance=1e-7)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    perm_strata <- rep(1, nrow(pheno2)) # avoid the stratified permutations
    names(perm_strata) <- rownames(pheno2)
    operm <- scan1perm(pr, pheno2, kinship_loco, addcovar=sex, Xcovar=Xcovar, n_perm=3,
                       perm_Xsp=TRUE, chr_lengths=chr_lengths(map),
                       perm_strata=perm_strata)
    expected <- list(A=cbind(liver= c(1.18430970191976991930, 0.85274796222960802528,  0.76981903233390269747),
                             spleen=c(1.97990260877063839470, 2.49362417238307543244,  1.25357633693210002157)),
                     X=cbind(liver= c(1.38737185599694012870, 0.40916780292085591642,  1.64260630582905675645, 0.73029802514640995703,
                                      1.48936109005232153457, 1.76554289216658433230, -0.19556097227382970849),
                             spleen=c(1.39377789099751669970, 1.24497337936097429711,  0.94154476906363604449, 0.50568976180393732101,
                                      0.70117350390854737974, 0.46582839192008501650,  0.98324516513467219436)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected, tolerance=1e-7)

})

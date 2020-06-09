context("scan1 permutations with binary phenotype")

iron <- read_cross2(system.file("extdata","iron.zip", package="qtl2"))
iron <- iron[,c(18,19,"X")]
map <- insert_pseudomarkers(iron$gmap, step=1)
pr <- calc_genoprob(iron, map, err=0.002)
kinship <- calc_kinship(pr)
kinship_loco <- calc_kinship(pr, "loco")

pheno1 <- apply(iron$pheno, 2, function(a) as.numeric(a > quantile(a, 0.7)))
rownames(pheno1) <- rownames(iron$pheno)
pheno2 <- pheno1; pheno2[5,1] <- NA

Xcovar <- get_x_covar(iron)
sex <- as.numeric(iron$covar$sex=="m")
names(sex) <- rownames(iron$covar)
perm_strata <- mat2strata(Xcovar)

test_that("scan1 permutations work (regression test)", {

    seed <- 3025685
    RNGkind("L'Ecuyer-CMRG")

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", n_perm=3)
    expected <- cbind(liver= c(1.9400869101173725539, 1.3826531402385029423, 2.32771714636574245105),
                      spleen=c(1.7406845630378597889, 1.0816176320165595826, 0.86215751014691477394))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    skip_on_cran()

    # no covariates or missing data
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", n_perm=3, perm_strata=perm_strata)
    expected <- cbind(liver= c(7.7800367429250911755, 8.2464978987470800575, 8.4332725697338304371),
                      spleen=c(2.6682023473596672147, 3.1865109293844255944, 3.9334182919580200632))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.0627090966702610331, 1.3019431062016479927, 1.4887177771883841615),
                      spleen=c(1.2907987947091044134, 1.0942134968225047942, 1.5908715555718089263))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.6988140471307104917, 3.3492519382653398452, 1.4887177771883841615),
                      spleen=c(1.3563064772388884194, 2.3361493956548571305, 1.5908715555718089263))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # sex and additive covariates + X-sp perms
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
    expected <- list(A=cbind(liver= c(1.0627090966702610331,  0.54091905189068256732, 0.51646113492044776194),
                             spleen=c(1.2907987947091044134,  1.09421349682250479418, 0.97567385090644620504)),
                     X=cbind(liver= c(0.92947674285197479094, 0.82005313512853206248, 0.41809257354587714417, 0.51654137988325032893,
                                      0.37354377722229514802, 1.8944354102049061339,  0.64614008414334023200),
                             spleen=c(1.86312726596119659916, 0.42851226191818625466, 1.82692521628932524891, 1.54255046971590559224,
                                      1.17331067477128669907, 1.9748233131450518840,  0.14820359491564261134)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno1, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3,
                        perm_Xsp=TRUE, chr_lengths=chr_lengths(map), perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, no covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, model="binary", n_perm=3)
    expected <- cbind(liver= c(1.9396621996084775219, 1.3860820684821959503, 2.39451070079611838537),
                      spleen=c(1.7406845630378597889, 1.0816176320165595826, 0.86215751014691477394))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3)
    expected <- cbind(liver= c(1.0607794866188555716, 1.2669252823191214929, 1.4888030146012454225),
                      spleen=c(1.2907987947091044134, 1.0942134968225047942, 1.5908715555718089263))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, plus sex interactive
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, model="binary", addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3)
    expected <- cbind(liver= c(1.6885642297968530556, 3.4770987617121420499, 1.4888030146012454225),
                      spleen=c(1.3563064772388884194, 2.3361493956548571305, 1.5908715555718089263))
    class(expected) <- c("scan1perm", "matrix")
    expect_equal(operm, expected)
    set.seed(seed)
    operm <- scan1perm(pr, pheno2, model="binary", addcovar=sex, Xcovar=Xcovar, intcovar=sex, n_perm=3, perm_strata=perm_strata)
    expect_equal(operm, expected)

    # one missing phenotype, sex and X-chr covariates, X-sp perms
    set.seed(seed)
    perm_strata <- rep(1, nrow(pheno2)) # avoid the stratified permutations
    names(perm_strata) <- rownames(pheno2)
    operm <- scan1perm(pr, pheno2, model="binary", addcovar=sex, Xcovar=Xcovar, n_perm=3,
                       perm_Xsp=TRUE, chr_lengths=chr_lengths(map),
                       perm_strata=perm_strata)
    expected <- list(A=cbind(liver= c(1.7660810516288876215, 0.40146200349977334554, 0.41346547186287807563),
                             spleen=c(1.5238369101838173947, 1.21467164163189522696, 0.79071777746629834382)),
                     X=cbind(liver= c(1.1856908571482307480, 1.6393912628303013435,  1.2664516688820839363, 0.72714128449122483744,
                                      2.7436059094899150068, 1.3434021611743531821,  0.49535669351789124448),
                             spleen=c(1.4926905975660105241, 1.2193200324768156406,  3.1386909220195633452, 0.56076286804034225497,
                                      1.4488076063209831545, 0.3586987681249098614,  1.20528954216223382900)))
    L <- c(A=61.200000000000002842,X=28.399999999999998579)
    attr(L, "is_x_chr") <- c(A=FALSE, X=TRUE)
    attr(expected, "chr_lengths") <- L
    class(expected) <- c("scan1perm", "list")
    expect_equal(operm, expected)

})

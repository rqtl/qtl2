context("I/O of yaml/json control files")

test_that("write_control_file gives correct output for riself", {

    dir <- tempdir()
    file <- file.path(dir, "grav2.yaml")
    on.exit(unlink(file))

    # write file
    write_control_file(file,
                    crosstype="riself",
                    geno_file="grav2_geno.csv",
                    gmap_file="grav2_gmap.csv",
                    pheno_file="grav2_pheno.csv",
                    phenocovar_file="grav2_phenocovar.csv",
                    geno_codes=c(L=1L, C=2L),
                    alleles=c("L", "C"),
                    na.strings=c("-", "NA"))

    # read it back in
    control <-  read_control_file(file)

    expected <- structure(list(crosstype = "riself", geno = "grav2_geno.csv",
                               pheno = "grav2_pheno.csv", phenocovar = "grav2_phenocovar.csv",
                               gmap = "grav2_gmap.csv", alleles = c("L", "C"),
                               genotypes = structure(list(L = 1L, C = 2L), .Names = c("L", "C")),
                               na.strings = c("-", "NA")),
                          .Names = c("crosstype", "geno", "pheno", "phenocovar",
                          "gmap", "alleles", "genotypes", "na.strings"))

    expect_equal(sort(names(control)), sort(names(expected)))
    expected <- expected[names(control)]
    expect_equal(control, expected)

})


test_that("write_control_file in JSON for riself", {

    dir <- tempdir()
    file <- file.path(dir, "grav2.json")
    on.exit(unlink(file))

    # write file
    write_control_file(file,
                    crosstype="riself",
                    geno_file="grav2_geno.csv",
                    gmap_file="grav2_gmap.csv",
                    pheno_file="grav2_pheno.csv",
                    phenocovar_file="grav2_phenocovar.csv",
                    geno_codes=c(L=1L, C=2L),
                    alleles=c("L", "C"),
                    na.strings=c("-", "NA"))

    # read it back in
    control <-  read_control_file(file)

    expected <- structure(list(crosstype = "riself", geno = "grav2_geno.csv",
                               pheno = "grav2_pheno.csv", phenocovar = "grav2_phenocovar.csv",
                               gmap = "grav2_gmap.csv", alleles = c("L", "C"),
                               genotypes = structure(list(L = 1L, C = 2L), .Names = c("L", "C")),
                               na.strings = c("-", "NA")),
                          .Names = c("crosstype", "geno", "pheno", "phenocovar",
                          "gmap", "alleles", "genotypes", "na.strings"))

    expect_equal(sort(names(control)), sort(names(expected)))
    expected <- expected[names(control)]
    expect_equal(control, expected)

})


test_that("write_control_file gives correct output for intercross", {

    dir <- tempdir()
    file <- file.path(dir, "iron.yaml")
    on.exit(unlink(file))

    # write file
    write_control_file(file,
                       crosstype="f2",
                       geno_file="iron_geno.csv",
                       gmap_file="iron_gmap.csv",
                       pheno_file="iron_pheno.csv",
                       covar_file="iron_covar.csv",
                       phenocovar_file="iron_phenocovar.csv",
                       geno_codes=c(SS=1L, SB=2L, BB=3L),
                       sex_covar="sex",
                       sex_codes=c(f="female", m="male"),
                       crossinfo_covar="cross_direction",
                       crossinfo_codes=c("(SxB)x(SxB)"=0L, "(BxS)x(BxS)"=1L),
                       xchr="X",
                       alleles=c("S", "B"),
                       na.strings=c("-", "NA"))

    # read it back in
    control <-  read_control_file(file)

    expected <- structure(list(crosstype = "f2", geno = "iron_geno.csv", pheno = "iron_pheno.csv",
                               phenocovar = "iron_phenocovar.csv", covar = "iron_covar.csv",
                               gmap = "iron_gmap.csv", alleles = c("S", "B"),
                               genotypes = structure(list(SS = 1L, SB = 2L, BB = 3L), .Names = c("SS", "SB", "BB")),
                               sex = structure(list(covar = "sex", f = "female", m = "male"),
                               .Names = c("covar", "f", "m")),
                               cross_info = structure(list(covar = "cross_direction", `(SxB)x(SxB)` = 0L, `(BxS)x(BxS)` = 1L),
                               .Names = c("covar", "(SxB)x(SxB)", "(BxS)x(BxS)")),
                               x_chr = "X", na.strings = c("-", "NA")),
                          .Names = c("crosstype", "geno", "pheno", "phenocovar",
                          "covar", "gmap", "alleles", "genotypes", "sex", "cross_info",
                          "x_chr", "na.strings"))

    expect_equal(sort(names(control)), sort(names(expected)))
    expected <- expected[names(control)]
    expect_equal(control, expected)

})


test_that("write_control_file in JSON intercross", {

    dir <- tempdir()
    file <- file.path(dir, "iron.json")
    on.exit(unlink(file))

    # write file
    write_control_file(file,
                       crosstype="f2",
                       geno_file="iron_geno.csv",
                       gmap_file="iron_gmap.csv",
                       pheno_file="iron_pheno.csv",
                       covar_file="iron_covar.csv",
                       phenocovar_file="iron_phenocovar.csv",
                       geno_codes=c(SS=1L, SB=2L, BB=3L),
                       sex_covar="sex",
                       sex_codes=c(f="female", m="male"),
                       crossinfo_covar="cross_direction",
                       crossinfo_codes=c("(SxB)x(SxB)"=0L, "(BxS)x(BxS)"=1L),
                       xchr="X",
                       alleles=c("S", "B"),
                       na.strings=c("-", "NA"))

    # read it back in
    control <-  read_control_file(file)

    expected <- structure(list(crosstype = "f2", geno = "iron_geno.csv", pheno = "iron_pheno.csv",
                               phenocovar = "iron_phenocovar.csv", covar = "iron_covar.csv",
                               gmap = "iron_gmap.csv", alleles = c("S", "B"),
                               genotypes = structure(list(SS = 1L, SB = 2L, BB = 3L), .Names = c("SS", "SB", "BB")),
                               sex = structure(list(covar = "sex", f = "female", m = "male"),
                               .Names = c("covar", "f", "m")),
                               cross_info = structure(list(covar = "cross_direction", `(SxB)x(SxB)` = 0L, `(BxS)x(BxS)` = 1L),
                               .Names = c("covar", "(SxB)x(SxB)", "(BxS)x(BxS)")),
                               x_chr = "X", na.strings = c("-", "NA")),
                          .Names = c("crosstype", "geno", "pheno", "phenocovar",
                          "covar", "gmap", "alleles", "genotypes", "sex", "cross_info",
                          "x_chr", "na.strings"))

    expect_equal(sort(names(control)), sort(names(expected)))
    expected <- expected[names(control)]
    expect_equal(control, expected)

})

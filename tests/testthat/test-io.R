context("input/output")

test_that("can read grav2 data", {

    zip_file <- system.file("extdata", "grav2.zip", package="qtl2geno")
    grav2 <- read_cross2(zip_file)

    # check that it contains the same stuff
    expect_equal(sort(names(grav2)),
                 c("alleles", "cross_info", "crosstype", "geno", "gmap", "is_female",
                   "is_x_chr", "pheno", "phenocovar"))

    # check summary
    expected <- structure(list(crosstype = "riself", nind = 162L, nind_geno=162L, nind_pheno=162L, nind_gnp=162L,
                               nchr = 5L, nmar = structure(c(26, 42, 64, 35, 67),
                                          .Names = c("1", "2", "3", "4", "5")),
                               npheno = 241L, ncovar = 0, nphenocovar = 1L,
                               totmar = 234), .Names = c("crosstype", "nind", "nind_geno", "nind_pheno", "nind_gnp",
                                              "nchr", "nmar", "npheno", "ncovar", "nphenocovar", "totmar"),
                          class = c("summary.cross2", "list"))
    expect_equal(summary(grav2), expected)

    # calculate QTL genotype probabilities
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    pr <- calc_genoprob(grav2, map)

})

test_that("can read iron data", {

    zip_file <- system.file("extdata", "iron.zip", package="qtl2geno")
    iron <- read_cross2(zip_file)

    # check that it contains the same stuff
    expect_equal(sort(names(iron)),
                 c("alleles", "covar", "cross_info", "crosstype", "geno", "gmap",
                   "is_female", "is_x_chr", "pheno", "phenocovar", "pmap"))

    # check summary
    expected <- structure(list(crosstype = "f2", nind = 284L, nind_geno = 284L, nind_pheno = 284L, nind_gnp = 284L,
                               nchr = 20L, nmar = structure(c(3, 5, 2, 2, 2, 2, 7, 8, 5,
                                           2, 7, 2, 2, 2, 2, 5, 2, 2, 2, 2),
                                           .Names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                           "11", "12", "13", "14", "15", "16", "17", "18", "19", "X")),
                               npheno = 2L, ncovar = 2L, nphenocovar = 1L, totmar = 66),
                          .Names = c("crosstype", "nind", "nind_geno", "nind_pheno", "nind_gnp", "nchr", "nmar", "npheno", "ncovar",
                          "nphenocovar", "totmar"), class = c("summary.cross2", "list"))
    expect_equal(summary(iron), expected)


    # calculate QTL genotype probabilities
    map <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, map)

})

test_that("read_pheno works", {

    # iron data
    ironfile <- system.file("extdata", "iron.zip", package="qtl2geno")
    dir <- tempdir()

    # read full data
    iron <- read_cross2(ironfile)

    # unzip data
    unzipped_files <- utils::unzip(ironfile, exdir=dir)

    # names of data files to be used
    phefile <- unzipped_files[grep("_pheno\\.csv$", unzipped_files)]
    phecovfile <- unzipped_files[grep("_phenocovar\\.csv$", unzipped_files)]

    # zip files to be created
    phezipfile <- paste0(phefile, ".zip")
    phecovzipfile <- paste0(phecovfile, ".zip")
    bothzipfile <- paste0(phefile, "_both.zip")

    # clean up
    on.exit(unlink(c(unzipped_files, phezipfile, phecovzipfile, bothzipfile)))

    # create zip files (-j: don't store file name, -q: be quiet)
    zip(phezipfile, phefile, flags="-j -q")
    zip(phecovzipfile, phecovfile, flags="-j -q")
    zip(bothzipfile, c(phefile, phecovfile), flags="-j -q")

    # read pheno as plain file
    expect_equal(read_pheno(phefile), iron$pheno)

    # read pheno as zip file
    expect_equal(read_pheno(phezipfile), iron$pheno)
    unzip(phezipfile, exdir=dirname(phezipfile)) # restore file

    # read pheno + phenocovar each in plain files
    phelist <- read_pheno(phefile, phecovfile)
    expect_equal(phelist$pheno, iron$pheno)
    expect_equal(phelist$phenocovar, iron$phenocovar)

    # read pheno + phenocovar in one zip file
    expect_error(read_pheno(bothzipfile))
    phelist <- read_pheno(bothzipfile, basename(phecovfile))
    expect_equal(phelist$pheno, iron$pheno)
    expect_equal(phelist$phenocovar, iron$phenocovar)
    unzip(bothzipfile, exdir=dirname(bothzipfile)) # restore files

    # read pheno from one zip file and phenocovar from plain file
    phelist <- read_pheno(phezipfile, phecovfile)
    expect_equal(phelist$pheno, iron$pheno)
    expect_equal(phelist$phenocovar, iron$phenocovar)
    unzip(phezipfile, exdir=dirname(phezipfile)) # restore file

    # read pheno from one zip file and phenocovar from another
    phelist <- read_pheno(phezipfile, phecovzipfile)
    expect_equal(phelist$pheno, iron$pheno)
    expect_equal(phelist$phenocovar, iron$phenocovar)
    unzip(phezipfile, exdir=dirname(phezipfile)) # restore file
    unzip(phecovzipfile, exdir=dirname(phecovzipfile)) # restore file

    # read pheno from zip file and phenocovar from plain file
    phelist <- read_pheno(phefile, phecovzipfile)
    expect_equal(phelist$pheno, iron$pheno)
    expect_equal(phelist$phenocovar, iron$phenocovar)
    unzip(phecovzipfile, exdir=dirname(phecovzipfile)) # restore file

})

test_that("Create zip file works", {

    # unzip iron data to a temporary directory
    ironfile <- system.file("extdata", "iron.zip", package="qtl2geno")
    dir <- tempdir()
    unzipped_files <- utils::unzip(ironfile, exdir=dir)
    on.exit(unlink(unzipped_files)) # clean up

    # zip the files
    zip_datafiles(file.path(dir, "iron.yaml"))
    on.exit(unlink(file.path(dir, "iron.zip")), add=TRUE)

    # file created?
    zipfile <- file.path(dir, "iron.zip")
    expect_true( file.exists(zipfile) )

    # unzip to tmp dir
    tmp_dir <- file.path(dir, "tmp")
    new_unzipped_files <- utils::unzip(zipfile, exdir=tmp_dir)
    on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE) # clean up

    # sample files as originally?
    ofiles <- sort( basename( unzipped_files ) )
    nfiles <- sort( basename( new_unzipped_files ) )
    expect_equal( ofiles, nfiles )

    # read them?
    x <- read_cross2(file.path(dir, "iron.yaml"))
    y <- read_cross2(file.path(dir, "iron.zip"))
    expect_equal(y, x)

})

context("zip data files")

test_that("zip_datafiles() works", {

    skip_if(isnt_karl(), "this test only run locally")

    # put files in temp directory
    tmpdir <- file.path(tempdir(), "test_zip_datafiles")

    # grav2 zip file, unzip into temp directory
    orig_zip_file <- system.file("extdata", "grav2.zip", package="qtl2")
    unzip(orig_zip_file, exdir=tmpdir)

    # read original version
    grav2a <- read_cross2(orig_zip_file)

    # create new version
    new_zip_file <- zip_datafiles(file.path(tmpdir, "grav2.yaml"), overwrite=TRUE)

    # load newly created file
    grav2b <- read_cross2(new_zip_file)
    expect_equal(grav2a, grav2b)

    # zip into another directory
    other_dir <- file.path(tempdir(), "test_zip_file_again")
    if(!dir.exists(other_dir)) {
        dir.create(other_dir)
    }
    zip_file <- file.path(other_dir, "other_zip.zip")
    zip_datafiles(file.path(tmpdir, "grav2.yaml"), zip_file, overwrite=TRUE)
    grav2c <- read_cross2(zip_file)
    expect_equal(grav2a, grav2c)

    # clean up
    unlink(other_dir, recursive=TRUE)
    unlink(tmpdir, recursive=TRUE)

})

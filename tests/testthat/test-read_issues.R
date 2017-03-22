context("read_cross2 handling of file inconsistencies")

test_that("read_cross2 deals with missing marker info", {

    ironfile <- system.file("extdata", "iron.zip", package="qtl2geno")
    iron <- read_cross2(ironfile)

    # unzip iron data to a temporary directory
    dir <- tempdir()
    unzipped_files <- utils::unzip(ironfile, exdir=dir)

    # drop two markers from the geno file
    set.seed(79646034)
    geno_file <- file.path(dir, "iron_geno.csv")
    g <- read_csv(geno_file, rownames_included=FALSE)
    dropped_index <- sort( sample(2:ncol(g), 2) )
    dropped_markers <- colnames(g)[dropped_index]
    write.table(g[,-dropped_index], file=geno_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    yaml_file <- file.path(dir, "iron.yaml")
    expect_warning( iron_sub <- read_cross2(yaml_file) )

    expect_equal(iron_sub, drop_markers(iron, dropped_markers))

    # drop two markers from the gmap file
    gmap_file <- file.path(dir, "iron_gmap.csv")
    gmap <- read_csv(gmap_file, rownames_included=FALSE)
    dropped_index2 <- sample((1:nrow(gmap))[-(dropped_index-1)], 2)
    dropped_markers2 <- gmap[dropped_index2,1]
    write.table(gmap[-dropped_index2,], file=gmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_sub2 <- read_cross2(yaml_file) ) # this gives an error at the moment

    expect_equal(iron_sub2, drop_markers(iron_sub, dropped_markers2))

    dropped_index <- c(dropped_index-1, dropped_index2)
    dropped_markers <- c(dropped_markers, dropped_markers2)

    # drop two markers from the pmap file
    pmap_file <- file.path(dir, "iron_pmap.csv")
    pmap <- read_csv(pmap_file, rownames_included=FALSE)
    dropped_index3 <- sample((1:nrow(gmap))[-dropped_index], 2)
    dropped_markers3 <- pmap[dropped_index3,1]
    write.table(pmap[-dropped_index3,], file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_sub3 <- read_cross2(yaml_file) ) # this gives an error at the moment

    expect_equal(iron_sub3, drop_markers(iron, c(dropped_markers, dropped_markers3)))


})

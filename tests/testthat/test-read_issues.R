context("read_cross2 handling of file inconsistencies")

test_that("read_cross2 deals with missing marker info", {

    ironfile <- system.file("extdata", "iron.zip", package="qtl2geno")
    iron <- read_cross2(ironfile)

    # unzip iron data to a temporary directory
    dir <- tempdir()
    unzipped_files <- utils::unzip(ironfile, exdir=dir)

    # drop two markers from the geno file
    geno_file <- file.path(dir, "iron_geno.csv")
    g <- read_csv(geno_file, rownames_included=FALSE)
    dropped_index <- c(20,28)
    dropped_markers <- colnames(g)[dropped_index]
    write.table(g[,-dropped_index], file=geno_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    yaml_file <- file.path(dir, "iron.yaml")
    expect_warning( iron_sub <- read_cross2(yaml_file) )

    expect_equal(iron_sub, drop_markers(iron, dropped_markers))

    # drop two markers from the gmap file
    # and give another two missing positions
    gmap_file <- file.path(dir, "iron_gmap.csv")
    gmap <- read_csv(gmap_file, rownames_included=FALSE)
    dropped_index2 <- c(14, 51)
    dropped_markers2 <- gmap[dropped_index2,1]
    na_index2 <- c(5, 56)
    na_markers2 <- gmap[na_index2,1]
    gmap[na_index2,3] <- NA
    write.table(gmap[-dropped_index2,], file=gmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_sub2 <- read_cross2(yaml_file) ) # this gives an error at the moment

    expect_equal(iron_sub2, drop_markers(iron_sub, c(dropped_markers2, na_markers2)))

    dropped_index <- c(dropped_index-1, dropped_index2, na_index2)
    dropped_markers <- c(dropped_markers, dropped_markers2, na_markers2)

    # drop two markers from the pmap file
    # and give another one a missing position
    pmap_file <- file.path(dir, "iron_pmap.csv")
    pmap <- read_csv(pmap_file, rownames_included=FALSE)
    dropped_index3 <- c(4, 29)
    dropped_markers3 <- pmap[dropped_index3,1]
    na_index3 <- 33
    na_markers3 <- pmap[na_index3,1]
    pmap[na_index3,3] <- NA
    write.table(pmap[-dropped_index3,], file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_sub3 <- read_cross2(yaml_file) ) # this gives an error at the moment
    expect_equal(iron_sub3, drop_markers(iron, c(dropped_markers, dropped_markers3, na_markers3)))

    # markers out of order in genotypes
    g <- read_csv(geno_file)
    g <- g[,sample(ncol(g)),drop=FALSE]
    write.table(cbind(marker=rownames(g), g), file=geno_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_randgeno <- read_cross2(yaml_file) )
    expect_equal(iron_sub3, iron_randgeno)

    # markers out of order in physical map
    pmap <- read_csv(pmap_file, rownames_included=FALSE)
    pmap <- pmap[sample(nrow(pmap)), , drop=FALSE]
    write.table(pmap, file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_randpmap <- read_cross2(yaml_file) )
    expect_equal(iron_sub3, iron_randpmap)



})

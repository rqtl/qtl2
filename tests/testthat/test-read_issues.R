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
    pmap_rand <- pmap
    pmap_rand <- pmap[sample(nrow(pmap)), , drop=FALSE]
    write.table(pmap_rand, file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_randpmap <- read_cross2(yaml_file) )
    expect_equal(iron_sub3, iron_randpmap)

    # markers out of order within a chromosome on genetic map
    gmap <- read_csv(gmap_file, rownames_included=FALSE)
    gmap_reordered <- gmap
    gmap_reordered[6:7,1] <- gmap[7:6,1] # swap two marker names
    write.table(gmap_reordered, file=gmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_marorder <- read_cross2(yaml_file) )
    expected <- iron_sub3
    expected$pmap[[2]] <- expected$pmap[[2]][c(2,1,3)]
    names(expected$gmap[[2]]) <- names(expected$gmap[[2]])[c(2,1,3)]
    expected$geno[[2]] <- expected$geno[[2]][,c(2,1,3)]
    expect_equal(iron_marorder, expected)

    # markers out of order within a chromosome on physical map
    ### first, put gmap and pmap back in order
    write.table(gmap, file=gmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(pmap, file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)
    pmap_reordered <- pmap
    neworder <- c(40,37,43,41,38,42,39)
    pmap_reordered[37:43,1] <- pmap[neworder,1] # reorder marker names on chr 11
    write.table(pmap_reordered, file=pmap_file, sep=",",
                row.names=FALSE, col.names=TRUE, quote=FALSE)

    expect_warning( iron_marorder <- read_cross2(yaml_file) )
    expected <- iron_sub3
    tmp <- expected$pmap[[11]]
    names(tmp) <- names(tmp)[c(4,1,7,5,2,6,3)]
    expected$pmap[[11]] <- tmp[names(expected$pmap[[11]])]
    expect_equal(iron_marorder, expected)

})

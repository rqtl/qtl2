context("check cross")

test_that("check_crosstype works appropriately", {

    good_types <- c("bc", "f2", "riself", "risib", "dh", "haploid")
    for(type in good_types)
        expect_true(check_crosstype(type))

    expect_error(check_crosstype("f2pk"))
    expect_false(check_crosstype("f2pk", FALSE))

    expect_error(check_crosstype("bada1gakadsf"))
    expect_false(check_crosstype("bada1gakadsf", FALSE))

})

test_that("count_invalid_genotypes works appropriately", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    expect_true(check_cross2(grav2))
    count <- count_invalid_genotypes(grav2)
    expect_true(all(count==0))

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    expect_true(check_cross2(iron))
    count <- count_invalid_genotypes(iron)
    expect_true(all(count==0))

    # add some errors
    set.seed(79648025)
    g <- grav2$geno[[1]]
    err <- sample(0:5, prod(dim(g)), replace=TRUE, prob=c(0.95, rep(0.01, 5)))
    grav2$geno[[1]] <- g+err
    count <- count_invalid_genotypes(grav2)
    expect_equivalent(count[,1] , rowSums(g+err > 2))
    expect_true(all(count[,-1] == 0))

    # add some errors to the other data
    g <- iron$geno[[1]]
    err <- sample(0:5, prod(dim(g)), replace=TRUE, prob=c(0.95, rep(0.01, 5)))
    iron$geno[[1]] <- g+err
    count <- count_invalid_genotypes(iron)
    expect_equivalent(count[,1] , rowSums(g+err > 5))
    expect_true(all(count[,-1] == 0))

})

test_that("check_cross2 gives proper warnings", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # missing geno
    iron_bad <- iron; iron_bad$geno <- NULL
    expect_warning(check_cross2(iron_bad))

    # missing gmap
    iron_bad <- iron; iron_bad$gmap <- NULL
    expect_warning(check_cross2(iron_bad))

    # missing cross_info
    iron_bad <- iron; iron_bad$cross_info <- NULL
    expect_warning(check_cross2(iron_bad))

    # missing is_female
    iron_bad <- iron; iron_bad$is_female <- NULL
    expect_warning(check_cross2(iron_bad))

    # missing is_x_chr
    iron_bad <- iron; iron_bad$is_x_chr <- NULL
    expect_warning(check_cross2(iron_bad))

    # missing alleles okay
    iron_bad <- iron; iron_bad$alleles <- NULL
    expect_true(check_cross2(iron_bad))

    # missing pmap okay
    iron_bad <- iron; iron_bad$pmap <- NULL
    expect_true(check_cross2(iron_bad))

    # gmap has no names()
    iron_bad <- iron; names(iron_bad$gmap) <- NULL
    expect_warning(check_cross2(iron_bad))

    # pmap has no names()
    iron_bad <- iron; names(iron_bad$pmap) <- NULL
    expect_warning(check_cross2(iron_bad))

    # gmap names != geno names
    iron_bad <- iron; names(iron_bad$gmap) <- 1:20
    expect_warning(check_cross2(iron_bad))

    # pmap names != geno names
    iron_bad <- iron; names(iron_bad$pmap) <- 1:20
    expect_warning(check_cross2(iron_bad))

})


test_that("check_cross2 works for MPP data", {

    if(isnt_karl()) skip("this test only run locally")

    do <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Recla/recla.zip")
    expect_true(check_cross2(do))

    # no founder genotypes
    do_bad <- do
    do_bad$founder_geno <- NULL
    expect_warning(check_cross2(do_bad))

    # no names in founder genotypes
    do_bad <- do
    names(do_bad$founder_geno) <- NULL
    expect_warning(check_cross2(do_bad))

    # mismatch in names of founder genotypes
    do_bad <- do
    names(do_bad$founder_geno) <- 1:20
    expect_warning(check_cross2(do_bad))

    # missing a founder in founder_geno
    do_bad <- do
    do_bad$founder_geno[[10]] <- do$founder_geno[[10]][-4,]
    expect_warning(check_cross2(do_bad))

    # missing a marker in founder_geno
    do_bad <- do
    do_bad$founder_geno[[12]] <- do$founder_geno[[12]][,-50]
    expect_warning(check_cross2(do_bad))

    # markers in founder_geno out of order
    do_bad <- do
    colnames(do_bad$founder_geno[[13]])[50:52] <- colnames(do$founder_geno[[13]])[52:50]
    expect_warning(check_cross2(do_bad))

    ### HS need founder geno
    hs <- do
    hs$crosstype <- "hs"
    expect_true(check_cross2(hs))

    hs_bad <- hs
    hs_bad$founder_geno <- NULL
    expect_warning(check_cross2(hs_bad))

    ### AIL3
    ail3 <- do
    ail3$crosstype <- "ail3"
    for(i in seq_along(ail3$founder_geno))
        ail3$founder_geno[[i]] <- ail3$founder_geno[[i]][1:3,,drop=FALSE]
    expect_true(check_cross2(ail3))

    ail3_bad <- ail3
    ail3_bad$founder_geno <- NULL
    expect_warning(check_cross2(ail3_bad))

    ### riself4
    riself4 <- do
    riself4$crosstype <- "riself4"
    for(x in c("geno", "founder_geno", "is_x_chr", "gmap", "pmap"))
        riself4[[x]] <- riself4[[x]][1:19]
    riself4$cross_info <- cbind(riself4$cross_info, 1L, 2L, 3L, 4L)[,-1,drop=FALSE]
    for(i in seq_along(riself4$founder_geno)) {
        riself4$founder_geno[[i]] <- riself4$founder_geno[[i]][1:4,,drop=FALSE]
        riself4$geno[[i]][riself4$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(riself4))

    riself4_bad <- riself4
    riself4_bad$founder_geno <- NULL
    expect_warning(check_cross2(riself4_bad))

    ### riself8
    riself8 <- do
    riself8$crosstype <- "riself8"
    for(x in c("geno", "founder_geno", "is_x_chr", "gmap", "pmap"))
        riself8[[x]] <- riself8[[x]][1:19]
    riself8$cross_info <- cbind(riself8$cross_info, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)[,-1,drop=FALSE]
    for(i in seq_along(riself8$founder_geno)) {
        riself8$geno[[i]][riself8$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(riself8))

    riself8_bad <- riself8
    riself8_bad$founder_geno <- NULL
    expect_warning(check_cross2(riself8_bad))

    ### riself16
    riself16 <- do
    riself16$crosstype <- "riself16"
    for(x in c("geno", "founder_geno", "is_x_chr", "gmap", "pmap"))
        riself16[[x]] <- riself16[[x]][1:19]
    riself16$cross_info <- cbind(riself16$cross_info, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L,
                                 9L, 10L, 11L, 12L, 13L, 14L, 15L, 16L)[,-1,drop=FALSE]
    for(i in seq_along(riself16$founder_geno)) {
        riself16$founder_geno[[i]] <- rbind(riself16$founder_geno[[i]], riself16$founder_geno[[i]])
        riself16$geno[[i]][riself16$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(riself16))

    riself16_bad <- riself16
    riself16_bad$founder_geno <- NULL
    expect_warning(check_cross2(riself16_bad))

    ### risib4
    risib4 <- do
    risib4$crosstype <- "risib4"
    risib4$cross_info <- cbind(risib4$cross_info, 1L, 2L, 3L, 4L)[,-1,drop=FALSE]
    for(i in seq_along(risib4$founder_geno)) {
        risib4$founder_geno[[i]] <- risib4$founder_geno[[i]][1:4,,drop=FALSE]
        risib4$geno[[i]][risib4$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(risib4))

    risib4_bad <- risib4
    risib4_bad$founder_geno <- NULL
    expect_warning(check_cross2(risib4_bad))

    ### risib8
    risib8 <- do
    risib8$crosstype <- "risib8"
    risib8$cross_info <- cbind(risib8$cross_info, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L)[,-1,drop=FALSE]
    for(i in seq_along(risib8$founder_geno)) {
        risib8$geno[[i]][risib8$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(risib8))

    risib8_bad <- risib8
    risib8_bad$founder_geno <- NULL
    expect_warning(check_cross2(risib8_bad))

    ### dh6
    dh6 <- do
    dh6$crosstype <- "dh6"
    for(x in c("geno", "founder_geno", "is_x_chr", "gmap", "pmap"))
        dh6[[x]] <- dh6[[x]][1:19]
    for(i in seq_along(dh6$founder_geno)) {
        dh6$founder_geno[[i]] <- dh6$founder_geno[[i]][1:6,]
        dh6$geno[[i]][dh6$geno[[i]] == 2] <- 0
    }
    expect_true(check_cross2(dh6))

    dh6_bad <- dh6
    dh6_bad$founder_geno <- NULL
    expect_warning(check_cross2(dh6_bad))

    ### dof1
    dof1 <- do
    dof1$crosstype <- "dof1"
    for(i in seq_along(dof1$founder_geno)) {
        dof1$founder_geno[[i]] <- dof1$founder_geno[[i]][c(1:8, 5),,drop=FALSE]
    }
    expect_true(check_cross2(dof1))

    dof1_bad <- dof1
    dof1_bad$founder_geno <- NULL
    expect_warning(check_cross2(dof1_bad))

})

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


test_that("check_cross2 works for DO data", {

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

})

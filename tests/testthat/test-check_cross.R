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

context("subset cross2, calc_genoprob, sim_geno")

test_that("subset.cross2 works (riself)", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    chr <- c(3,5)
    grav2sub <- grav2[,chr]

    expected <- grav2
    expected$geno <- expected$geno[chr]
    expected$gmap <- expected$gmap[chr]
    expected$is_x_chr <- expected$is_x_chr[chr]

    expect_equal(grav2sub, expected)

    ind <- c("79", "100", "123", "160")
    grav2subA <- grav2sub[ind,]
    grav2subB <- grav2[ind,chr]

    for(i in seq(along=expected$geno))
        expected$geno[[i]] <- expected$geno[[i]][ind,,drop=FALSE]
    expected$pheno <- expected$pheno[ind,,drop=FALSE]
    expected$cross_info <- expected$cross_info[ind,,drop=FALSE]
    expected$is_female <- expected$is_female[ind]

    expect_equal(grav2subA, expected)
    expect_equal(grav2subB, expected)

})

test_that("subset.cross2 works (F2)", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    chr <- rep(FALSE, 20)
    chr[c(4, 20)] <- TRUE
    ironsub <- iron[,chr]

    expected <- iron
    expected$geno <- expected$geno[chr]
    expected$gmap <- expected$gmap[chr]
    expected$is_x_chr <- expected$is_x_chr[chr]

    expect_equal(ironsub, expected)

    ind <- c("79", "100", "123", "160", "284")
    ironsubA <- ironsub[ind,]
    ironsubB <- iron[ind,chr]

    for(i in seq(along=expected$geno))
        expected$geno[[i]] <- expected$geno[[i]][ind,,drop=FALSE]
    expected$pheno <- expected$pheno[ind,,drop=FALSE]
    expected$covar <- expected$covar[ind,,drop=FALSE]
    expected$cross_info <- expected$cross_info[ind,,drop=FALSE]
    expected$is_female <- expected$is_female[ind]

    expect_equal(ironsubA, expected)
    expect_equal(ironsubB, expected)

})

test_that("subset.calc_genoprob works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(4,"X")]

    pr <- calc_genoprob(ironsub, step=5, err=0.01)
    prsub <- pr[,"X"]

    at <- attributes(pr)
    expected <- unclass(pr)["X"]
    attr(expected, "is_x_chr") <- at$is_x_chr["X"]
    attr(expected, "map") <- at$map["X"]
    attr(expected, "crosstype") <- at$crosstype
    attr(expected, "cross_info") <- at$cross_info
    attr(expected, "class") <- at$class
    attr(expected, "alleles") <- at$alleles

    expect_equal(prsub, expected)

    ind <- c("5", "50", "55", "280")
    prsub <- pr[ind, "X"]

    expected[["X"]] <- expected[["X"]][ind,,,drop=FALSE]
    attr(expected, "cross_info") <- at$cross_info[ind,,drop=FALSE]

    expect_equal(prsub, expected)

})

test_that("subset.sim_geno works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(4,"X")]

    dr <- sim_geno(iron, step=5, err=0.01)
    drsub <- dr[,"X"]

    at <- attributes(dr)
    expected <- unclass(dr)["X"]
    attr(expected, "is_x_chr") <- at$is_x_chr["X"]
    attr(expected, "map") <- at$map["X"]
    attr(expected, "map") <- at$map["X"]
    attr(expected, "crosstype") <- at$crosstype
    attr(expected, "cross_info") <- at$cross_info
    attr(expected, "class") <- at$class

    expect_equal(drsub, expected)

    ind <- c("5", "50", "55", "280")
    drsub <- dr[ind, "X"]

    expected[["X"]] <- expected[["X"]][ind,,,drop=FALSE]
    attr(expected, "cross_info") <- at$cross_info[ind,,drop=FALSE]

    expect_equal(drsub, expected)

})

test_that("subset.calc_genoprob works with reduction to grid and/or allele prob", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(5,8,"X")]

    pr <- calc_genoprob(ironsub, step=2.5, err=0.002)

    pr_grid <- probs_to_grid(pr)
    pr_a <- genoprob_to_alleleprob(pr)

    set.seed(20150918)
    ind <- sample(rownames(iron$geno[[1]]), 10, replace=FALSE)
    chr <- c("5", "X")
    pr_sub <- pr[ind,chr]

    pr_grid_sub <- pr_grid[ind,chr]
    pr_sub_grid <- probs_to_grid(pr_sub)
    expect_equal(pr_grid_sub, pr_sub_grid)

    pr_a_sub <- pr_a[ind,chr]
    pr_sub_a <- genoprob_to_alleleprob(pr_sub)
    expect_equal(pr_a_sub, pr_sub_a)

})

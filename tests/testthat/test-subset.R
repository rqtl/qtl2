context("subset cross2, calc_genoprob, sim_geno, viterbi")

test_that("subset.cross2 works (riself)", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))

    # possible problems in indexes
    expect_error(grav2[,c(18:19,"X")])
    expect_warning(grav2[,c(1,2,18)])
    expect_equal(grav2[,1:2], suppressWarnings(grav2[,c(1:2,18)]))
    expect_error(grav2[201:250,])
    expect_warning(grav2[101:250,])
    expect_equal(grav2[101:162,], suppressWarnings(grav2[101:250,]))

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
    expected$pmap <- expected$pmap[chr]
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

    # test negative subscripts
    expect_equal(iron[,"-X"], iron[,1:19])
    expect_equal(iron[,"-4"], iron[,c(1:3,5:19,"X")])
    expect_equal(iron[c("-189", "-190"),], iron[-c(189,190),])

})

test_that("subset.calc_genoprob works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(4,"X")]

    map <- insert_pseudomarkers(iron$gmap, step=5)

    pr <- calc_genoprob(iron, map, err=0.01)
    prsub <- pr[,"X"]

    expected <- unclass(pr)["X"]
    class(expected) <- class(pr)
    attr(expected, "is_x_chr") <- attr(pr, "is_x_chr")["X"]
    for(obj in c("alleles", "alleleprobs", "crosstype"))
        attr(expected, obj) <- attr(pr, obj)

    expect_equal(prsub, expected)

    ind <- c("5", "50", "55", "280")
    prsub <- pr[ind, "X"]

    expected[["X"]] <- expected[["X"]][ind,,,drop=FALSE]

    expect_equal(prsub, expected)

})

test_that("subset.sim_geno works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(4,"X")]

    map <- insert_pseudomarkers(iron$gmap, step=5)

    dr <- sim_geno(iron, map, err=0.01)
    drsub <- dr[,"X"]

    expected <- unclass(dr)["X"]
    class(expected) <- class(dr)
    attr(expected, "is_x_chr") <- attr(dr, "is_x_chr")["X"]
    for(obj in c("alleles", "crosstype"))
        attr(expected, obj) <- attr(dr, obj)

    expect_equal(drsub, expected)

    ind <- c("5", "50", "55", "280")
    drsub <- dr[ind, "X"]

    expected[["X"]] <- expected[["X"]][ind,,,drop=FALSE]

    expect_equal(drsub, expected)

})

test_that("subset.calc_genoprob works with reduction to grid and/or allele prob", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(5,8,"X")]

    map <- insert_pseudomarkers(ironsub$gmap, step=2.5)

    pr <- calc_genoprob(ironsub, map, err=0.002)

    grid <- calc_grid(ironsub$gmap, step=2.5)
    pr_grid <- probs_to_grid(pr, grid)
    pr_a <- genoprob_to_alleleprob(pr)

    set.seed(20150918)
    ind <- sample(rownames(iron$geno[[1]]), 10, replace=FALSE)
    chr <- c("5", "X")
    pr_sub <- pr[ind,chr]

    pr_grid_sub <- pr_grid[ind,chr]
    pr_sub_grid <- probs_to_grid(pr_sub, grid[chr])
    expect_equal(pr_grid_sub, pr_sub_grid)

    pr_a_sub <- pr_a[ind,chr]
    pr_sub_a <- genoprob_to_alleleprob(pr_sub)
    expect_equal(pr_a_sub, pr_sub_a)

})

test_that("subset.viterbi works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    ironsub <- iron[,c(4,"X")]

    map <- insert_pseudomarkers(iron$gmap, step=5)
    g <- viterbi(iron, map, err=0.01)
    gsub <- g[,"X"]

    expected <- unclass(g)["X"]
    class(expected) <- class(g)
    attr(expected, "is_x_chr") <- attr(g, "is_x_chr")["X"]
    attr(expected, "crosstype") <- attr(g, "crosstype")
    attr(expected, "alleles") <- attr(g, "alleles")

    expect_equal(gsub, expected)

    ind <- c("5", "50", "55", "280")
    gsub <- g[ind, "X"]

    expected[["X"]] <- expected[["X"]][ind,,drop=FALSE]

    expect_equal(gsub, expected)

})

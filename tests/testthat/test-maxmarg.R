context("maxmarg")

test_that("maxmarg works for F2 data", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,8] # only chr 8

    pr <- calc_genoprob(iron, step=2.5, err=0.01)

    z <- maxmarg(pr, minprob=0)
    g <- z$geno[[1]]
    expected <- apply(pr$probs[[1]], c(1,3), function(a) seq(along=a)[a==max(a)])
    expect_equal(g, expected)

    # single locus
    zz <- maxmarg(pr, minprob=0, chr="8", pos=45.5)
    expect_equal(zz, g[,"D8Mit40"])

    # repeat with minprob=0.95
    z <- maxmarg(pr, minprob=0.95)
    g <- z$geno[[1]]
    expected <- apply(pr$probs[[1]], c(1,3), function(a) { b <- max(a); ifelse(b > 0.95, which.max(a), NA) })
    expect_equal(g, expected)

    # single locus
    zz <- maxmarg(pr, minprob=0.95, chr="8", pos=45.5)
    expect_equal(zz, g[,"D8Mit40"])

    # repeat with return_value="genotypes"
    zchar <- maxmarg(pr, minprob=0.95, return_char=TRUE)
    gchar <- g;gchar[,1:ncol(g)] <- c("SS", "SB", "BB")[g]
    expect_equal(zchar$geno[[1]], gchar)

    # single locus
    zzchar <- maxmarg(pr, minprob=0.95, chr="8", pos=45.5, return_char=TRUE)
    gg <- c("SS", "SB", "BB")[zz]
    names(gg) <- names(zz)
    expect_equal(zzchar, gg)
})

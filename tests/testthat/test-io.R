
context("input/output")

test_that("can read grav2 data", {

    zip_file <- system.file("extdata", "grav2.zip", package="qtl2")
    grav2 <- read_cross2(zip_file)

    # attempt to calculate QTL genotype probabilities on chr 1
    pmap <- insert_pseudomarkers(grav2$gmap[[1]], step=1, stepwidth="max", pmar_stem="c1.loc")
    rf <- mf.h(diff(pmap))
    pr_chr1 <- calc_genoprob(grav2$crosstype, t(grav2$geno[[1]]), grav2$is_x_chr[1],
                             grav2$is_female, t(grav2$cross_info),
                             mf.h(diff(pmap)), attr(pmap, "index")-1, 0.01)

})

test_that("can read iron data", {

    zip_file <- system.file("extdata", "iron.zip", package="qtl2")
    iron <- read_cross2(zip_file)

    # attempt to calculate QTL genotype probabilities on X chr
    pmap <- insert_pseudomarkers(iron$gmap[["X"]], step=1, pmar_stem="cX.loc")
    rf <- mf.h(diff(pmap))
    pr_chrX <- calc_genoprob(iron$crosstype, t(iron$geno[["X"]]), iron$is_x_chr["X"],
                             iron$is_female, t(iron$cross_info),
                             mf.h(diff(pmap)), attr(pmap, "index")-1, 0.01)

})

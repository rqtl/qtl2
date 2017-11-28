context("is_same for comparing objects")

test_that("is_same works for calc_genoprob results", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2[6:12,], map, error_prob=0.002)

    expect_true(is_same(probs, probs))

    for(i in seq(along=probs))
        expect_true(is_same(probs[[i]], probs[[i]]))

    expect_false(is_same(probs, probs[1:2,]))
})

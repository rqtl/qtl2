context("Scale kinship matrix")

test_that("scale_kinship works for RIL", {

    set.seed(49265251)

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    sim <- calc_kinship(probs)

    scaled <- scale_kinship(sim)

    # row and colnames okay
    expect_equal(dimnames(scaled), dimnames(sim))

    # check some values
    expect_equal(diag(scaled), setNames(rep(1, nrow(sim)), colnames(sim)))

    d <- sqrt(diag(sim))
    n <- 10
    pairs <- matrix(sample(ncol(sim), 2*n, replace=TRUE), ncol=2)
    for(i in 1:n)
        expect_equal(scaled[pairs[i,1], pairs[i,2]],
                     setNames(sim[pairs[i,1], pairs[i,2]] / (d[pairs[i,1]] * d[pairs[i,2]]), NULL) )


})

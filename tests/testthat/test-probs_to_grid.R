context("Reduce probabilities to pseudomarker grid")

test_that("probs_to_grid works", {

    # try it out
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    orig_dim <- sapply(probs, dim)
    grid <- calc_grid(grav2$gmap, step=1)
    probs_sub <- probs_to_grid(probs, grid)
    new_dim <- sapply(probs_sub, dim)

    # reduced dimension match what we'd expect?
    npmar <- 1 + floor(sapply(map, function(a) diff(range(a))))
    expected_dim <- orig_dim; expected_dim[3,] <- npmar
    expect_equal(new_dim, expected_dim)

    # test results
    expected <- probs
    for(i in seq(along=probs)) {
        expected[[i]] <- probs[[i]][,,grid[[i]],drop=FALSE]
    }

    expect_equal(probs_sub, expected)

})

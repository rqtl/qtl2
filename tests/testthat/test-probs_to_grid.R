context("Reduce probabilities to pseudomarker grid")

test_that("probs_to_grid works", {

    # try it out
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    orig_dim <- sapply(probs$probs, dim)
    probs_sub <- probs_to_grid(probs)
    new_dim <- sapply(probs_sub$probs, dim)

    # map at which calculations were done
    map <- probs$map

    # reduced dimension match what we'd expect?
    npmar <- 1 + floor(sapply(map, function(a) diff(range(a))))
    expected_dim <- orig_dim; expected_dim[3,] <- npmar
    expect_equal(new_dim, expected_dim)

    # test results
    expected <- probs
    for(i in seq(along=probs$chrID)) {
        grid <- probs$grid[[i]]
        expected$probs[[i]] <- probs$probs[[i]][,,grid,drop=FALSE]

        map[[i]] <- map[[i]][grid]
    }

    expected$map <- map
    expected$grid <- NULL
    expect_equal(probs_sub, expected)
})

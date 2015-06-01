context("Reduce probabilities to pseudomarker grid")

test_that("probs_to_grid works", {

    # try it out
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    orig_dim <- sapply(probs, dim)
    probs_sub <- probs_to_grid(probs)
    new_dim <- sapply(probs_sub, dim)

    # map at which calculations were done
    map <- attr(probs, "map")

    # reduced dimension match what we'd expect?
    npmar <- 1 + floor(sapply(map, function(a) diff(range(a))))
    expected_dim <- orig_dim; expected_dim[2,] <- npmar
    expect_equal(new_dim, expected_dim)

    # test results
    expected <- probs
    for(i in seq(along=probs)) {
        mapat <- attributes(map[[i]])
        grid <- mapat$grid
        expected[[i]] <- probs[[i]][,grid,,drop=FALSE]

        map[[i]] <- map[[i]][grid]
        for(j in c("grid", "index"))
            mapat[[j]] <- mapat[[j]][grid]
        for(j in names(mapat)[names(mapat) != "names"])
            attr(map[[i]], j) <- mapat[[j]]
    }

    attr(expected, "map") <- map
    expect_equal(probs_sub, expected)
})

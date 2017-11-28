context("interpolation of genotype probabilities")

test_that("interp_genoprob works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[1:20,c("1", "2", "X")]
    probs <- calc_genoprob(iron, iron$gmap, error_prob=0.002)

    # no added positions
    expect_equal(interp_genoprob(probs, iron$gmap), probs)

    # mess up the order
    map <- iron$gmap
    map[[2]] <- map[[2]][c(1,3,2,4,5)]
    map[[2]][2:3] <- c(48.1, 56.8)
    expect_error(interp_genoprob(probs, map))

    # add a couple of pseudomarkers
    map <- insert_pseudomarkers(iron$gmap, step=20, stepwidth="max")
    expected <- probs
    expected[[1]] <- array(dim=c(nrow(probs[[1]]), ncol(probs[[1]]), length(map[[1]])))
    dimnames(expected[[1]]) <- c(dimnames(probs[[1]])[1:2], list(names(map[[1]])))
    expected[[1]][,,c(1,3,6)] <- probs[[1]]
    expected[[1]][,,2] <- (probs[[1]][,,1] + probs[[1]][,,2])/2
    p <- (map[[1]][4]-map[[1]][3])/(map[[1]][6]-map[[1]][3])
    expected[[1]][,,4] <- probs[[1]][,,2]*(1-p) + probs[[1]][,,3]*p
    p <- (map[[1]][5]-map[[1]][3])/(map[[1]][6]-map[[1]][3])
    expected[[1]][,,5] <- probs[[1]][,,2]*(1-p) + probs[[1]][,,3]*p
    expected[[3]] <- array(dim=c(nrow(probs[[3]]), ncol(probs[[3]]), length(map[[3]])))
    dimnames(expected[[3]]) <- c(dimnames(probs[[3]])[1:2], list(names(map[[3]])))
    expected[[3]][,,c(1,3)] <- probs[[3]]
    expected[[3]][,,2] <- (probs[[3]][,,1] + probs[[3]][,,2])/2

    expect_equal(interp_genoprob(probs, map), expected)

    # add some positions off the ends
    map[[1]] <- c(off_c1=10, map[[1]])
    map[[2]] <- c(off_c2a=10, off_c2b=20, map[[2]], off_c2c=80)
    map[[3]] <- c(off_cXa=10, map[[3]], off_cXb=60)

    expected2 <- expected
    for(i in 1:3) {
        expected2[[i]] <- array(dim=c(nrow(probs[[i]]), ncol(probs[[i]]), length(map[[i]])))
        dimnames(expected2[[i]]) <- c(dimnames(probs[[i]])[1:2], list(names(map[[i]])))
    }
    expected2[[1]][,,1:7] <- expected[[1]][,,c(1,1:6)]
    expected2[[2]][,,1:8] <- expected[[2]][,,c(1,1,1:5,5)]
    expected2[[3]][,,1:5] <- expected[[3]][,,c(1,1:3,3)]

    expect_equal(interp_genoprob(probs, map), expected2)

})

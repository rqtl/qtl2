context("Conversion of genotype probabilities to allele probabilities")

test_that("genoprob_to_alleleprob works for RIL", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    probs <- calc_genoprob(grav2, step=1, error_prob=0.002)
    allele_probs <- genoprob_to_alleleprob(probs)

    attr(probs, "alleles") <- TRUE # should be only difference between the two
    expect_equal(allele_probs, probs)

})

test_that("genoprob_to_alleleprob works for F2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    probs <- calc_genoprob(iron, step=1, error_prob=0.002)
    allele_probs <- genoprob_to_alleleprob(probs)

    expected <- probs
    for(i in which(!attr(probs, "is_x_chr"))) { # loop over autosomes
        expected[[i]] <- probs[[i]][,1:2,]
        expected[[i]][,1,] <- probs[[i]][,1,] + probs[[i]][,2,]/2
        expected[[i]][,2,] <- probs[[i]][,3,] + probs[[i]][,2,]/2
    }
    attr(expected, "alleles") <- TRUE

    expect_equal(allele_probs, expected)

})

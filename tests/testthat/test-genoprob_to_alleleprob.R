context("Conversion of genotype probabilities to allele probabilities")

test_that("genoprob_to_alleleprob works for RIL", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    allele_probs <- genoprob_to_alleleprob(probs)

    # expected result, hardly changed
    expected <- probs
    attr(expected, "alleleprobs") <- TRUE
    for(i in seq(along=probs))
        colnames(expected[[i]]) <- c("L", "C")

    expect_equal(allele_probs, expected)

})

test_that("genoprob_to_alleleprob works for F2", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)
    allele_probs <- genoprob_to_alleleprob(probs)

    f2_geno2alle <-
        function(prob, x_chr=FALSE)
        {
            if(x_chr) {
                cn <- colnames(prob)
                prob[,1,] <- prob[,1,]+prob[,2,]/2+prob[,3,]/2 + prob[,5,]
                prob[,2,] <- prob[,4,]+prob[,2,]/2+prob[,3,]/2 + prob[,6,]
                prob <- prob[,1:2,]
                colnames(prob) <- substr(cn[c(1,3)], 1, 1)
                return(prob)

            } else {
                cn <- colnames(prob)
                prob[,1,] <- prob[,1,]+prob[,2,]/2
                prob[,2,] <- prob[,3,]+prob[,2,]/2
                prob <- prob[,1:2,]
                colnames(prob) <- substr(cn[c(1,3)], 1, 1)
                return(prob)
            }
        }

    expected <- probs
    is_x_chr <- attr(probs, "is_x_chr")
    for(i in seq(along=probs)) # loop over chromosomes
        expected[[i]] <- f2_geno2alle(probs[[i]], is_x_chr[i])
    attr(expected, "alleleprobs") <- TRUE

    expect_equal(allele_probs, expected)

})

test_that("genoprob_to_alleleprob works when multi-core", {
    if(isnt_karl()) skip("this test only run locally")

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(iron$gmap, step=1)
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    allele_probs <- genoprob_to_alleleprob(probs)
    allele_probs_mc <- genoprob_to_alleleprob(probs, cores=4)
    expect_equal(allele_probs_mc, allele_probs)

    # following shouldn't really matter, since no calculations are done
    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    map <- insert_pseudomarkers(grav2$gmap, step=1)
    probs <- calc_genoprob(grav2, map, error_prob=0.002)
    allele_probs <- genoprob_to_alleleprob(probs)
    allele_probs_mc <- genoprob_to_alleleprob(probs, cores=4)
    expect_equal(allele_probs_mc, allele_probs)

})

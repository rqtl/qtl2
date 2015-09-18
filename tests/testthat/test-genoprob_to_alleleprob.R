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

    f2_geno2alle <-
        function(prob, x_chr=FALSE)
        {
            if(x_chr) {
                prob[,1,] <- prob[,1,]+prob[,2,]/2+prob[,3,]/2
                prob[,2,] <- prob[,4,]+prob[,2,]/2+prob[,3,]/2
                return(prob[,-(3:4),])

            } else {
                prob[,1,] <- prob[,1,]+prob[,2,]/2
                prob[,2,] <- prob[,3,]+prob[,2,]/2
                return(prob[,1:2,])
            }
        }

    expected <- probs
    is_x_chr <- attr(probs, "is_x_chr")
    for(i in seq(along=probs)) # loop over chromosomes
        expected[[i]] <- f2_geno2alle(probs[[i]], is_x_chr[i])
    attr(expected, "alleles") <- TRUE

    expect_equal(allele_probs, expected)

})

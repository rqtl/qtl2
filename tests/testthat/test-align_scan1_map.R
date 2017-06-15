context("align_scan1_map")

test_that("align_scan1_map works", {

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))

    # insert pseudomarkers into map
    map <- insert_pseudomarkers(iron$gmap, step=1)

    # calculate genotype probabilities
    probs <- calc_genoprob(iron, map, error_prob=0.002)

    # grab phenotypes and covariates; ensure that covariates have names attribute
    pheno <- iron$pheno[,1]
    covar <- match(iron$covar$sex, c("f", "m")) # make numeric
    names(covar) <- rownames(iron$covar)

    Xcovar <- get_x_covar(iron)

    # genome scan
    out <- scan1(probs, iron$pheno, Xcovar=Xcovar, addcovar=covar)

    # calculate coefficients for chromosome 7
    coef <- scan1coef(probs[,7], pheno, addcovar=covar)

    # not subseted
    result <- align_scan1_map(out, map)
    expect_equal(result$scan1, out)
    expect_equal(result$map, map)

    # chr 7 for map
    result <- align_scan1_map(out, map[7])
    expect_equal( result$scan1, subset(out, map, chr="7") )
    expect_equal( result$map, map[7] )

    # chr 7 for scan
    result <- align_scan1_map(subset(out,map,chr="7"), map)
    expect_equal( result$scan1, subset(out, map, chr="7") )
    expected_map <- map[7]
    attr(expected_map, "is_x_chr") <- attr(map, "is_x_chr")[7]
    expect_equal( result$map, expected_map )

    # coefficients
    result <- align_scan1_map(coef, map)
    expect_equal( result$scan1, coef )
    expected_map <- map[7]
    attr(expected_map, "is_x_chr") <- attr(map, "is_x_chr")[7]
    expect_equal( result$map, expected_map )

})

context("basic HMM functions in general RIL")

test_that("genril n_gen, n_alleles work", {

    expect_equal(nalleles("genril38"), 38)
    expect_equal(nalleles("genril6"), 6)

    expect_equal(test_ngen("genril38", FALSE), 38)
    expect_equal(test_ngen("genril8", FALSE), 8)
    expect_equal(test_ngen("genril38", TRUE), 38)
    expect_equal(test_ngen("genril8", TRUE), 8)

    # throw error for "genril1"
    expect_error(nalleles("genril1"))

})

test_that("genril possible_gen work", {

    expect_equal(test_possible_gen("genril38", FALSE, FALSE, 1:38), 1:38)

})

# FIX_ME: test check_geno

test_that("genril init work", {

    set.seed(20181105)
    alpha <- sample(1:10, 38, replace=TRUE)
    init<- alpha/sum(alpha)
    expect_equal( sapply(1:38, function(i) test_init("genril38", i, FALSE, FALSE, c(8,alpha))), log(init))

    expect_equal( sapply(1:38, function(i) test_init("genril38", i, FALSE, FALSE, c(8,alpha))), log(init))

})

# FIX_ME: test emit

test_that("genril step works", {

    nf <- 38
    alpha_int <- sample(1:10, nf, replace=TRUE)
    alpha <- alpha_int/sum(alpha_int)

    for(rf in c(0.01, 0.1, 0.45)) {
        for(ngen in c(3, 5)) {
            expected <- matrix(rep(alpha*(1-(1-rf)^ngen),nf), ncol=nf, byrow=TRUE)
            diag(expected) <- alpha + (1-alpha)*(1-rf)^ngen

            result <- matrix(ncol=nf, nrow=nf)
            for(i in 1:nf) {
                for(j in 1:nf) {
                    result[i,j] <- test_step(paste0("genril",nf), i, j, rf, FALSE, FALSE, c(ngen, alpha_int))
                }
            }

            expect_equal(result, log(expected))
        }
    }

})


test_that("genril geno_names work", {

    # if 38 founders, using upper case and lower case letters. Ugh.
    alleles <- c(LETTERS, letters[1:(38-26)])
    expect_equal( geno_names("genril38", alleles, FALSE), paste0(alleles, alleles) )

    # could also use two-letter allele codes, but ugly
    alleles <- c(paste0("A",LETTERS), paste0("B", LETTERS[1:(38-26)]))
    expect_equal( geno_names("genril38", alleles, FALSE), paste0(alleles, alleles) )

})

test_that("genril nrec work", {

    x <- matrix(ncol=38, nrow=38)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=38)

    res38 <- matrix(ncol=38, nrow=38)
    for(i in 1:38) {
        for(j in 1:38) {
            res38[i,j] <- test_nrec("genril38", i, j, FALSE, FALSE, c(3, rep(1, 38)))
        }
    }

    expect_equal( res38, x )

})

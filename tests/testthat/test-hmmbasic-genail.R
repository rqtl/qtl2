context("basic HMM functions in general AIL")

### TODO none of this has been revised yet

test_that("genail n_gen, n_alleles work", {

    expect_equal(nalleles("genail38"), 38)
    expect_equal(nalleles("genail6"), 6)

    expect_equal(test_ngen("genail38", FALSE), 38)
    expect_equal(test_ngen("genail8", FALSE), 8)
    expect_equal(test_ngen("genail38", TRUE), 38)
    expect_equal(test_ngen("genail8", TRUE), 8)

    # throw error for "genail1"
    expect_error(nalleles("genail1"))

})

test_that("genail possible_gen work", {

    expect_equal(test_possible_gen("genail38", FALSE, FALSE, 1:38), 1:38)

})

# FIX_ME: test check_geno

test_that("genail init work", {

    set.seed(20181105)
    alpha <- sample(1:10, 38, replace=TRUE)
    init<- alpha/sum(alpha)
    expect_equal( sapply(1:38, function(i) test_init("genail38", i, FALSE, FALSE, c(8,alpha))), log(init))

    expect_equal( sapply(1:38, function(i) test_init("genail38", i, FALSE, FALSE, c(8,alpha))), log(init))

})

# FIX_ME: test emit

test_that("genail step works", {

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
                    result[i,j] <- test_step(paste0("genail",nf), i, j, rf, FALSE, FALSE, c(ngen, alpha_int))
                }
            }

            cat(rf, ngen, "\n")
            expect_equal(result, log(expected))
        }
    }

})


test_that("genail geno_names work", {

    # if 38 founders, using upper case and lower case letters. Ugh.
    alleles <- c(LETTERS, letters[1:(38-26)])
    expect_equal( geno_names("genail38", alleles, FALSE), paste0(alleles, alleles) )

    # could also use two-letter allele codes, but ugly
    alleles <- c(paste0("A",LETTERS), paste0("B", LETTERS[1:(38-26)]))
    expect_equal( geno_names("genail38", alleles, FALSE), paste0(alleles, alleles) )

})

test_that("genail nrec work", {

    x <- matrix(ncol=38, nrow=38)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=38)

    res38 <- matrix(ncol=38, nrow=38)
    for(i in 1:38) {
        for(j in 1:38) {
            res38[i,j] <- test_nrec("genail38", i, j, FALSE, FALSE, c(3, rep(1, 38)))
        }
    }

    expect_equal( res38, x )

})

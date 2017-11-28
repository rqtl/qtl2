context("basic HMM functions in 6-way doubled haploids")

test_that("dh6 n_gen, n_alleles work", {

    expect_equal(nalleles("dh6"), 6)

    expect_equal(test_ngen("dh6", FALSE), 6)

})

test_that("dh6 possible_gen work", {

    expect_equal(test_possible_gen("dh6", FALSE, FALSE, 1:6), 1:6)

})

# FIX_ME: test check_geno

test_that("dh6 init work", {

    expect_equal( sapply(1:6, function(i) test_init("dh6", i, FALSE, FALSE, 1:6)), log(rep(1/6, 6)) )

})

# FIX_ME: test emit

test_that("dh6 step works", {

    for(rf in c(0.01, 0.1, 0.45)) {
        for(ngen in c(3, 5)) {
            expected <- matrix((5 - (5-6*rf)*(1-rf)^(ngen-2))/30, ncol=6, nrow=6)
            diag(expected) <- (1 + (5-6*rf)*(1-rf)^(ngen-2))/6

            result <- matrix(ncol=6, nrow=6)
            for(i in 1:6) {
                for(j in 1:6) {
                    result[i,j] <- test_step("dh6", i, j, rf, FALSE, FALSE, ngen)
                }
            }

            expect_equal(result, log(expected))
        }
    }

})


test_that("dh6 geno_names work", {

    expect_equal( geno_names("dh6", LETTERS[5:10], FALSE), paste0(LETTERS[5:10], LETTERS[5:10]) )

})

test_that("dh6 nrec work", {

    x <- matrix(ncol=6, nrow=6)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=6)

    res6 <- matrix(ncol=6, nrow=6)
    for(i in 1:6) {
        for(j in 1:6) {
            res6[i,j] <- test_nrec("dh6", i, j, FALSE, FALSE, 3)
        }
    }

    expect_equal( res6, x )

})

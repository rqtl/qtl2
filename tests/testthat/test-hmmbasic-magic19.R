context("basic HMM functions for magic19 (19-way RIL by selfing)")

test_that("magic19 n_gen, n_alleles work", {

    expect_equal(nalleles("magic19"), 19)

    expect_equal(test_ngen("magic19", FALSE), 19)

})

test_that("magic19 possible_gen work", {

    expect_equal(test_possible_gen("magic19", FALSE, FALSE, 1:19), 1:19)

})

# FIX_ME: test check_geno

test_that("magic19 init work", {

    expect_equal( sapply(1:19, function(i) test_init("magic19", i, FALSE, FALSE, integer(0))), log(rep(1/19, 19)) )

})

# FIX_ME: test emit

test_that("magic19 step works", {

    for(r in c(0.01, 0.1, 0.45)) {
        expected <- matrix(r*(90 - 54*r + 18*r^2)/(1+2*r)/19/18, ncol=19, nrow=19)
        diag(expected) <- (19 - 52*r + 54*r^2 - 18*r^3)/(1+2*r)/19

        result <- matrix(ncol=19, nrow=19)
        for(i in 1:19) {
            for(j in 1:19) {
                result[i,j] <- test_step("magic19", i, j, r, FALSE, FALSE, integer(0))
            }
        }

        expect_equal(result, log(expected))
    }

})



test_that("magic19 geno_names work", {

    expect_equal( geno_names("magic19", LETTERS[1:19], FALSE), paste0(LETTERS[1:19], LETTERS[1:19]) )

})

test_that("magic19 nrec work", {

    x <- matrix(ncol=19, nrow=19)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=19)

    res19 <- matrix(ncol=19, nrow=19)
    for(i in 1:19) {
        for(j in 1:19) {
            res19[i,j] <- test_nrec("magic19", i, j, FALSE, FALSE, integer(0))
        }
    }

    expect_equal( res19, x)

})

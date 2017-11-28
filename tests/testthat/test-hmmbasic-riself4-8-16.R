context("basic HMM functions in 4-, 8-, and 16-way RIL by selfing")

test_that("riself4-8-16 n_gen, n_alleles work", {

    expect_equal(nalleles("riself4"), 4)
    expect_equal(nalleles("riself8"), 8)
    expect_equal(nalleles("riself16"), 16)

    expect_equal(test_ngen("riself4", FALSE), 4)
    expect_equal(test_ngen("riself8", FALSE), 8)
    expect_equal(test_ngen("riself16", FALSE), 16)

})

test_that("riself4-8-16 possible_gen work", {

    expect_equal(test_possible_gen("riself4", FALSE, FALSE, 1:4), 1:4)
    expect_equal(test_possible_gen("riself8", FALSE, FALSE, 1:8), 1:8)
    expect_equal(test_possible_gen("riself16", FALSE, FALSE, 1:16), 1:16)

})

# FIX_ME: test check_geno

test_that("riself4-8-16 init work", {

    expect_equal( sapply(1:4, function(i) test_init("riself4", i, FALSE, FALSE, 1:4)), log(rep(1/4, 4)) )
    expect_equal( sapply(1:8, function(i) test_init("riself8", i, FALSE, FALSE, 1:8)), log(rep(1/8, 8)) )
    expect_equal( sapply(1:16, function(i) test_init("riself16", i, FALSE, FALSE, 1:16)), log(rep(1/16, 16)) )

})

# FIX_ME: test emit

test_that("riself4 step works", {

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/4/(1+2*rf), ncol=4, nrow=4)
        diag(expected) <- (1-rf)/4/(1+2*rf)

        result <- matrix(ncol=4, nrow=4)
        for(i in 1:4) {
            for(j in 1:4) {
                result[i,j] <- test_step("riself4", i, j, rf, FALSE, FALSE, 1:4)
            }
        }

        expect_equal(result, log(expected))
    }

})


test_that("riself8 step works", {

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/16/(1+2*rf), ncol=8, nrow=8)
        diag(expected) <- (1-rf)^2/8/(1+2*rf)
        block <- list(c(1,2), c(3,4), c(5,6), c(7,8))
        for(i in seq_along(block)) {
            b <- block[[i]]
            expected[b[1],b[2]] <- expected[b[2],b[1]] <- rf*(1-rf)/8/(1+2*rf)
        }

        result <- matrix(ncol=8, nrow=8)
        for(i in 1:8) {
            for(j in 1:8) {
                result[i,j] <- test_step("riself8", i, j, rf, FALSE, FALSE, 1:8)
            }
        }

        expect_equal(result, log(expected))
    }

    # test with different order of founders
    forder <- c(3, 8, 5, 2, 1, 7, 6, 4)
    forder_rev <- c(5,4,1,8,3,7,6,2)

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/16/(1+2*rf), ncol=8, nrow=8)
        diag(expected) <- (1-rf)^2/8/(1+2*rf)
        block <- list(c(1,2), c(3,4), c(5,6), c(7,8))
        for(i in seq_along(block)) {
            b <- block[[i]]
            expected[b[1],b[2]] <- expected[b[2],b[1]] <- rf*(1-rf)/8/(1+2*rf)
        }
        expected <- expected[forder_rev, forder_rev]

        result <- matrix(ncol=8, nrow=8)
        for(i in 1:8) {
            for(j in 1:8) {
                result[i,j] <- test_step("riself8", i, j, rf, FALSE, FALSE, forder)
            }
        }

        expect_equal(result, log(expected))
    }

})




test_that("riself16 step works", {

    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/64/(1+2*rf), ncol=16, nrow=16)
        expected[1:4,1:4] <- expected[5:8,5:8] <-
            expected[9:12,9:12] <- expected[13:16,13:16] <- rf*(1-rf)/32/(1+2*rf)
        diag(expected) <- (1-rf)^3/16/(1+2*rf)
        block <- list(c(1,2), c(3,4), c(5,6), c(7,8),c(9,10),c(11,12),c(13,14),c(15,16))
        for(i in seq_along(block)) {
            b <- block[[i]]
            expected[b[1],b[2]] <- expected[b[2],b[1]] <- rf*(1-rf)^2/16/(1+2*rf)
        }

        result <- matrix(ncol=16, nrow=16)
        for(i in 1:16) {
            for(j in 1:16) {
                result[i,j] <- test_step("riself16", i, j, rf, FALSE, FALSE, 1:16)
            }
        }

        expect_equal(result, log(expected))
    }

    # test with founders in a different order
    forder <- c(12, 5, 15, 13, 14, 1, 2, 3, 9, 10, 6, 11, 8, 4, 16, 7)
    forder_rev <- c(6,7,8,14,2,11,16,13,9,10,12,1,4,5,3,15)
    for(rf in c(0.01, 0.1, 0.45)) {
        expected <- matrix(rf/64/(1+2*rf), ncol=16, nrow=16)
        expected[1:4,1:4] <- expected[5:8,5:8] <-
            expected[9:12,9:12] <- expected[13:16,13:16] <- rf*(1-rf)/32/(1+2*rf)
        diag(expected) <- (1-rf)^3/16/(1+2*rf)
        block <- list(c(1,2), c(3,4), c(5,6), c(7,8),c(9,10),c(11,12),c(13,14),c(15,16))
        for(i in seq_along(block)) {
            b <- block[[i]]
            expected[b[1],b[2]] <- expected[b[2],b[1]] <- rf*(1-rf)^2/16/(1+2*rf)
        }
        expected <- expected[forder_rev, forder_rev]

        result <- matrix(ncol=16, nrow=16)
        for(i in 1:16) {
            for(j in 1:16) {
                result[i,j] <- test_step("riself16", i, j, rf, FALSE, FALSE, forder)
            }
        }

        expect_equal(result, log(expected))
    }

})




test_that("riself4-8-16 geno_names work", {

    expect_equal( geno_names("riself4", LETTERS[5:8], FALSE), paste0(LETTERS[5:8], LETTERS[5:8]) )
    expect_equal( geno_names("riself8", LETTERS[2:9], FALSE), paste0(LETTERS[2:9], LETTERS[2:9]) )
    expect_equal( geno_names("riself16", LETTERS[11:26], FALSE), paste0(LETTERS[11:26], LETTERS[11:26]) )

})

test_that("riself4-8-16 nrec work", {

    x <- matrix(ncol=16, nrow=16)
    x <- matrix(as.numeric(col(x) != row(x)), ncol=16)

    res4 <- matrix(ncol=4, nrow=4)
    res8 <- matrix(ncol=8, nrow=8)
    res16 <- matrix(ncol=16, nrow=16)
    for(i in 1:16) {
        for(j in 1:16) {
            if(i<5 && j<5) res4[i,j] <- test_nrec("riself4", i, j, FALSE, FALSE, 1:4)
            if(i<9 && j<9) res8[i,j] <- test_nrec("riself8", i, j, FALSE, FALSE, 1:8)
            res16[i,j] <- test_nrec("riself16", i, j, FALSE, FALSE, 1:16)
        }
    }

    expect_equal( res4, x[1:4,1:4])
    expect_equal( res8, x[1:8,1:8])
    expect_equal( res16, x)

})

context("basic HMM functions in phase-known heterogeneous stock")

test_that("Phase-known HS nalleles works", {
    expect_equal(nalleles("hspk"), 8)
})

test_that("Phase-known HS check_geno works", {

    # observed genotypes
    for(i in 0:5) {
        # Autosome
        expect_true(test_check_geno("hspk", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("hspk", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_true(test_check_geno("hspk", i, TRUE, TRUE, FALSE, 20))
    }
    for(i in c(-1, 6)) {
        # Autosome
        expect_false(test_check_geno("hspk", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("hspk", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_false(test_check_geno("hspk", i, TRUE, TRUE, FALSE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 1:64) {
        # Autosome
        expect_true(test_check_geno("hspk", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("hspk", i, FALSE, TRUE, TRUE, 20))
    }
    for(i in c(0, 65)) {
        # Autosome
        expect_false(test_check_geno("hspk", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("hspk", i, FALSE, TRUE, TRUE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 64 + 1:8) {
        expect_true(test_check_geno("hspk", i, FALSE, TRUE, FALSE, 20))
    }
    for(i in c(0, 1, 2, 63, 64, 64+9)) {
        expect_false(test_check_geno("hspk", i, FALSE, TRUE, FALSE, 20))
    }

})

test_that("Phase-known HS n_gen works", {

    expect_equal(test_ngen("hspk", FALSE), 64)
    expect_equal(test_ngen("hspk", TRUE),  64+8)

})

test_that("Phase-known HS possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("hspk", FALSE, FALSE, 20), 1:64)

    # X female
    expect_equal(test_possible_gen("hspk", TRUE, TRUE, 20), 1:64)

    # X male
    expect_equal(test_possible_gen("hspk", TRUE, FALSE, 20), 64+(1:8))

})

test_that("Phase-known HS init works", {

    # autosome and female X
    for(i in 1:36) {
        expect_equal(test_init("hspk", i, FALSE, FALSE, 20), log(1/64))
        expect_equal(test_init("hspk", i, TRUE,  TRUE,  20), log(1/64))
    }

    # male X
    for(i in 64 + 1:8)
        expect_equal(test_init("hspk", i, TRUE,  FALSE, 20), log(1/8))

})


test_that("Phase-known HS emit works", {

    fgen <- c(1,3,0,1,3,0,1,3) # founder genotypes; 0=missing, 1=AA, 3=BB
    err <- 0.01

    # Autosome or female X
    # truth = homA: AA (1), AD (7), AG (22), DD (10), DG (25), GG (28)
    expected <- log(c(1-err, err/2, err/2, 1-err/2, err))
    for(trueg in c(1,7,22,10,25,28)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = het: AB (2), AE (11), AH (29), BD (8), BG (23), DE (14), DH (32), EG (26), GH (35)
    expected <- log(c(err/2, 1-err, err/2, 1-err/2, 1-err/2))
    for(trueg in c(2,11,29,8,23,14,32,26,35)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = homB: BB (3), BE (12), BH (30), EE (15), EH (33), HH (36)
    expected <- log(c(err/2, err/2, 1-err, err, 1-err/2))
    for(trueg in c(3,12,30,15,33,36)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = A-: AC (4), AF (16), CD (9), DF (19), CG (24), FG (27)
    expected <- log(c(1-err,1,err,1-err,err))
    for(trueg in c(4,16,9,19,24,27)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = B-: BC (5), BF (17), CE (13), EF (20), CH (31), FH (34)
    expected <- log(c(err,1,1-err,err,1-err))
    for(trueg in c(5,17,13,20,31,34)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }

    # male X: treat het as missing
    # truth = hemA: A (1+64), D (4+64), G (7+64)
    expected <- log(c(1-err, 1, err, 1-err, err))
    for(trueg in 64+c(1,4,7))
        for(obsg in 1:5)
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = hemB: BB (2+64), E (5+64), H (8+64)
    expected <- log(c(err, 1, 1-err, err, 1-err))
    for(trueg in 64+c(2,5,8))
        for(obsg in 1:5)
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = missing: C (3+64), F (6+64)
    expected <- rep(0,5)
    for(trueg in 64+c(3,6))
        for(obsg in 1:5)
            expect_equal(test_emit("hspk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])


})

test_that("Phase-known HS step works", {

    # autosome
    ng <- 64
    trmat <- matrix(nrow=ng, ncol=ng)
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("hspk", gl, gr, rf, FALSE, FALSE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))
            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

    # female X
    ng <- 64
    trmat <- matrix(nrow=ng, ncol=ng)
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("hspk", gl, gr, rf, FALSE, FALSE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))
            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

    # male X
    ng <- 8
    trmat <- matrix(nrow=ng, ncol=ng)
    for(rf in c(0.01, 0.001, 0.0001)) {
        for(ngen in c(6, 12, 50)) {

            for(gl in 1:ng)
                for(gr in 1:ng)
                    trmat[gl,gr] <- test_step("do", 64+gl, 64+gr, rf, TRUE, FALSE, ngen)

            # no missing values
            expect_true(all(!is.na(trmat)))
            # all in (-Inf, 0]
            expect_true(all(trmat > -Inf & trmat <= 0))
            # rows sum to 1
            expect_equal( rowSums(exp(trmat)), rep(1, ng))
            # maximum value on the diagonal
            expect_equal( apply(trmat, 1, which.max), 1:ng)
        }
    }

})

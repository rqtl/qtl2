context("basic HMM functions in 3-way AIL")

test_that("AIL3 nalleles works", {
    expect_equal(nalleles("ail3"), 3)
})

test_that("AIL3 check_geno works", {

    # observed genotypes
    for(i in 0:5) {
        # Autosome
        expect_true(test_check_geno("ail3", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("ail3", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_true(test_check_geno("ail3", i, TRUE, TRUE, FALSE, 20))
    }
    for(i in c(-1, 6)) {
        # Autosome
        expect_false(test_check_geno("ail3", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("ail3", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_false(test_check_geno("ail3", i, TRUE, TRUE, FALSE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 1:6) {
        # Autosome
        expect_true(test_check_geno("ail3", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("ail3", i, FALSE, TRUE, TRUE, 20))
    }
    for(i in c(0, 7)) {
        # Autosome
        expect_false(test_check_geno("ail3", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("ail3", i, FALSE, TRUE, TRUE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 6 + 1:3) {
        expect_true(test_check_geno("ail3", i, FALSE, TRUE, FALSE, 20))
    }
    for(i in c(0:6, 6+4)) {
        expect_false(test_check_geno("ail3", i, FALSE, TRUE, FALSE, 20))
    }

})

test_that("AIL3 n_gen works", {

    expect_equal(test_ngen("ail3", FALSE), 6)
    expect_equal(test_ngen("ail3", TRUE),  9)

})

test_that("AIL3 possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("ail3", FALSE, FALSE, 20), 1:6)

    # X female
    expect_equal(test_possible_gen("ail3", TRUE, TRUE, 20), 1:6)

    # X male
    expect_equal(test_possible_gen("ail3", TRUE, FALSE, 20), 6+(1:3))

})

test_that("AIL3 init works", {

    hom <- cumsum(1:3)
    het <- (1:6)[!((1:6) %in% cumsum(1:3))]
    male <- 6 + (1:3)

    for(i in hom) {
        # autosome and female X
        expect_equal(test_init("ail3", i, FALSE, FALSE, 20), log(1/9))
        expect_equal(test_init("ail3", i, TRUE,  TRUE,  20), log(1/9))
    }

    for(i in het) {
        # autosome and female X
        expect_equal(test_init("ail3", i, FALSE, FALSE, 20), log(2/9))
        expect_equal(test_init("ail3", i, TRUE,  TRUE,  20), log(2/9))
    }

    for(i in male)
        expect_equal(test_init("ail3", i, TRUE,  FALSE, 20), log(1/3))

})


test_that("AIL3 emit works", {

    fgen <- c(1,3,0) # founder genotypes; 0=missing, 1=AA, 3=BB
    err <- 0.01

    # Autosome or female X
    # truth = homA: AA (1)
    expected <- log(c(1-err, err/2, err/2, 1-err/2, err))
    for(trueg in 1) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = het: AB (2)
    expected <- log(c(err/2, 1-err, err/2, 1-err/2, 1-err/2))
    for(trueg in 2) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = homB: BB (3)
    expected <- log(c(err/2, err/2, 1-err, err, 1-err/2))
    for(trueg in 3) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = A-: AC (4)
    expected <- log(c(1-err,1,err,1-err,err))
    for(trueg in 4) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = B-: BC (5)
    expected <- log(c(err,1,1-err,err,1-err))
    for(trueg in 5) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }

    # male X: treat het as missing
    # truth = hemA: A (1+6)
    expected <- log(c(1-err, 1, err, 1-err, err))
    for(trueg in 7)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = hemB: BB (2+6)
    expected <- log(c(err, 1, 1-err, err, 1-err))
    for(trueg in 8)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = missing: C (3+6)
    expected <- rep(0,5)
    for(trueg in 9)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])

})

test_that("AIL3 step works for autosome", {

    pAA <- function(n, r, k=3)
    {
        stopifnot(n >= 2)
        (1 - (1 - k + k*r)*(1-r)^(n-2))/k^2
    }

    rf <- 0.02
    ngen <- 10

    pAA <- pAA(ngen, rf)
    R <- (1 - 3*pAA)

    # expected values
    expected <- matrix(ncol=6, nrow=6)
    # AA->AA
    expected[1,1] <- expected[3,3] <- expected[6,6] <- (1-R)^2
    # AB->AB
    expected[2,2] <- expected[4,4] <- expected[5,5] <- (1-R)^2 + (R/2)^2
    # AA->BB
    expected[1,3] <- expected[1,6] <- expected[3,6] <- (R/2)^2
    # AA->AB
    expected[1,2] <- expected[1,4] <- expected[3,5] <- expected[2,3] <-
        expected[4,6] <- expected[5,6] <- (1-R)*(R/2)
    # AA->BC
    expected[1,5] <- expected[3,4] <- expected[2,6] <- (R/2)^2
    # AB->AC
    expected[2,4] <- expected[2,5] <- expected[4,5] <- (1-R)*(R/2) + (R/2)^2
    expected[lower.tri(expected)] <- t(expected)[lower.tri(expected)]
    expected <- log(expected)

    result <- matrix(ncol=6, nrow=6)
    for(i in 1:6)
        for(j in 1:6)
            result[i,j] <- test_step("ail3", i, j, rf, FALSE, FALSE, ngen)

    expect_equal(result, expected)

})


test_that("AIL3 step works for X chromosome", {

    pAA <- function(n, r, k=3)
    {
        stopifnot(n >= 2)
        z <- sqrt((1-r)*(9-r))

        male <- (1-r)/k * ( (1-r+z)/(2*z) * ((1-r-z)/4)^(n-2) +
                            (-1+r+z)/(2*z) * ((1-r+z)/4)^(n-2)) +
            (2-r)/(2*k) * ( (1-r-z)/2 * (1-r+z)/(2*z) * ((1-r-z)/4)^(n-2) +
                            (1-r+z)/2 * (-1+r+z)/(2*z) * ((1-r+z)/4)^(n-2)) +
            ( (r^2 + r*(z-5))/(k^2*(3+r+z)) * (1-r+z)/(2*z) * ((1-r-z)/4)^(n-2) +
              (r^2 - r*(z+5))/(k^2*(3+r-z)) * (-1+r+z)/(2*z) * ((1-r+z)/4)^(n-2) + 1/k^2)

        female <- (1-r)/k * ( (-1/z)*((1-r-z)/4)^(n-2) +
                              (1/z)*((1-r+z)/4)^(n-2)) +
            (2-r)/(2*k) * ( (1-r-z)/2 * (-1/z)*((1-r-z)/4)^(n-2) +
                            (1-r+z)/2 *  (1/z)*((1-r+z)/4)^(n-2)) +
            ( (r^2 + r*(z-5))/(k^2*(3+r+z)) * (-1/z)*((1-r-z)/4)^(n-2) +
              (r^2 - r*(z+5))/(k^2*(3+r-z)) * (1/z)*((1-r+z)/4)^(n-2)  + 1/k^2)

        if(any(r==0))
            male[r==0] <- female[r==0] <- 1/k

        if(any(r==1)) {
            if(n==2) male[r==1] <- 0
            if(n>2) male[r==1] <- 1/k^2
            if(n==2) female[r==1] <- 1/(2*k)
            if(n==3) female[r==1] <- 1/(2*k^2)
            if(n>3)  female[r==1] <- 1/k^2
        }

        list(female=female, male=male)
    }

    rf <- 0.02
    ngen <- 10

    pAA <- pAA(ngen, rf)

    # female X
    R <- (1 - 3*pAA$female)

    # expected values
    expected <- matrix(ncol=6, nrow=6)
    # AA->AA
    expected[1,1] <- expected[3,3] <- expected[6,6] <- (1-R)^2
    # AB->AB
    expected[2,2] <- expected[4,4] <- expected[5,5] <- (1-R)^2 + (R/2)^2
    # AA->BB
    expected[1,3] <- expected[1,6] <- expected[3,6] <- (R/2)^2
    # AA->AB
    expected[1,2] <- expected[1,4] <- expected[3,5] <- expected[2,3] <-
        expected[4,6] <- expected[5,6] <- (1-R)*(R/2)
    # AA->BC
    expected[1,5] <- expected[3,4] <- expected[2,6] <- (R/2)^2
    # AB->AC
    expected[2,4] <- expected[2,5] <- expected[4,5] <- (1-R)*(R/2) + (R/2)^2
    expected[lower.tri(expected)] <- t(expected)[lower.tri(expected)]
    expected <- log(expected)

    result <- matrix(ncol=6, nrow=6)
    for(i in 1:6)
        for(j in 1:6)
            result[i,j] <- test_step("ail3", i, j, rf, TRUE, TRUE, ngen)

    expect_equal(result, expected)

    # male X
    R <- (1 - 3*pAA$male)

    # expected values
    expected <- matrix(R/2, ncol=3, nrow=3)
    diag(expected) <- 1-R
    expected <- log(expected)

    result <- matrix(ncol=3, nrow=3)
    for(i in 1:3)
        for(j in 1:3)
            result[i,j] <- test_step("ail3", 6+i, 6+j, rf, TRUE, FALSE, ngen)

    expect_equal(result, expected)

})



test_that("geno_names works", {
    auto <- c("AA", "AB", "BB", "AC", "BC", "CC")
    X <- c(auto, paste0(LETTERS[1:3], "Y"))

    expect_equal(geno_names("ail3", LETTERS[1:3], FALSE), auto)
    expect_equal(geno_names("ail3", LETTERS[1:3], TRUE), X)
})

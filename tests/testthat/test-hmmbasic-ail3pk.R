context("basic HMM functions in phase-known 3-way AIL")

test_that("AIL3PK nalleles works", {
    expect_equal(nalleles("ail3pk"), 3)
})

test_that("AIL3PK check_geno works", {

    # observed genotypes
    for(i in 0:5) {
        # Autosome
        expect_true(test_check_geno("ail3pk", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("ail3pk", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_true(test_check_geno("ail3pk", i, TRUE, TRUE, FALSE, 20))
    }
    for(i in c(-1, 6)) {
        # Autosome
        expect_false(test_check_geno("ail3pk", i, TRUE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("ail3pk", i, TRUE, TRUE, TRUE, 20))
        # Male X
        expect_false(test_check_geno("ail3pk", i, TRUE, TRUE, FALSE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 1:9) {
        # Autosome
        expect_true(test_check_geno("ail3pk", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_true(test_check_geno("ail3pk", i, FALSE, TRUE, TRUE, 20))
    }
    for(i in c(0, 10)) {
        # Autosome
        expect_false(test_check_geno("ail3pk", i, FALSE, FALSE, FALSE, 20))
        # Female X
        expect_false(test_check_geno("ail3pk", i, FALSE, TRUE, TRUE, 20))
    }

    # true genotypes, autosome and female X
    for(i in 9 + 1:3) {
        expect_true(test_check_geno("ail3pk", i, FALSE, TRUE, FALSE, 20))
    }
    for(i in c(0:9, 9+4)) {
        expect_false(test_check_geno("ail3pk", i, FALSE, TRUE, FALSE, 20))
    }

})

test_that("AIL3PK n_gen works", {

    expect_equal(test_ngen("ail3pk", FALSE), 9)
    expect_equal(test_ngen("ail3pk", TRUE),  12)

})

test_that("AIL3PK possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("ail3pk", FALSE, FALSE, 20), 1:9)

    # X female
    expect_equal(test_possible_gen("ail3pk", TRUE, TRUE, 20), 1:9)

    # X male
    expect_equal(test_possible_gen("ail3pk", TRUE, FALSE, 20), 9+(1:3))

})

test_that("AIL3PK init works", {

    for(i in 1:9) {
        # autosome and female X
        expect_equal(test_init("ail3pk", i, FALSE, FALSE, 20), log(1/9))
        expect_equal(test_init("ail3pk", i, TRUE,  TRUE,  20), log(1/9))
    }

    for(i in 10:12)
        expect_equal(test_init("ail3pk", i, TRUE,  FALSE, 20), log(1/3))

})


test_that("AIL3PK emit works", {

    fgen <- c(1,3,0) # founder genotypes; 0=missing, 1=AA, 3=BB
    err <- 0.01

    # Autosome or female X
    # truth = homA: AA (1)
    expected <- log(c(1-err, err/2, err/2, 1-err/2, err))
    for(trueg in 1) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = het: AB (2) or BA (7)
    expected <- log(c(err/2, 1-err, err/2, 1-err/2, 1-err/2))
    for(trueg in c(2,7)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = homB: BB (3)
    expected <- log(c(err/2, err/2, 1-err, err, 1-err/2))
    for(trueg in 3) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = A-: AC (4) or CA (8)
    expected <- log(c(1-err,1,err,1-err,err))
    for(trueg in c(4,8)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }
    # truth = B-: BC (5) or CB (9)
    expected <- log(c(err,1,1-err,err,1-err))
    for(trueg in c(5,9)) {
        for(obsg in 1:5) {
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen, FALSE, FALSE, 20), expected[obsg])
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE,  TRUE, 20), expected[obsg])
        }
    }

    # male X: treat het as missing
    # truth = hemA: A (1+9)
    expected <- log(c(1-err, 1, err, 1-err, err))
    for(trueg in 10)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = hemB: BB (2+9)
    expected <- log(c(err, 1, 1-err, err, 1-err))
    for(trueg in 11)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])
    # truth = missing: C (3+9)
    expected <- rep(0,5)
    for(trueg in 12)
        for(obsg in 1:5)
            expect_equal(test_emit("ail3pk", obsg, trueg, err, fgen,  TRUE, FALSE, 20), expected[obsg])

})

test_that("AIL3PK step works for autosome", {

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
    expected <- matrix(1, ncol=9, nrow=9)
    a1 <- c("A", "A", "B", "A", "B", "C", "B", "C", "C")
    a2 <- c("A", "B", "B", "C", "C", "C", "A", "A", "B")

    for(i in 1:9) {
        for(j in 1:9) {
            expected[i,j] <- expected[i,j] * ifelse(a1[i]==a1[j], 1-R, R/2)
            expected[i,j] <- expected[i,j] * ifelse(a2[i]==a2[j], 1-R, R/2)
        }
    }
    expected <- log(expected)

    result <- matrix(ncol=9, nrow=9)
    for(i in 1:9)
        for(j in 1:9)
            result[i,j] <- test_step("ail3pk", i, j, rf, FALSE, FALSE, ngen)

    expect_equal(result, expected)

})


test_that("AIL3PK step works for X chromosome", {

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
    expected <- matrix(1, ncol=9, nrow=9)
    a1 <- c("A", "A", "B", "A", "B", "C", "B", "C", "C")
    a2 <- c("A", "B", "B", "C", "C", "C", "A", "A", "B")

    for(i in 1:9) {
        for(j in 1:9) {
            expected[i,j] <- expected[i,j] * ifelse(a1[i]==a1[j], 1-R, R/2)
            expected[i,j] <- expected[i,j] * ifelse(a2[i]==a2[j], 1-R, R/2)
        }
    }
    expected <- log(expected)

    result <- matrix(ncol=9, nrow=9)
    for(i in 1:9)
        for(j in 1:9)
            result[i,j] <- test_step("ail3pk", i, j, rf, TRUE, TRUE, ngen)

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
            result[i,j] <- test_step("ail3pk", 9+i, 9+j, rf, TRUE, FALSE, ngen)

    expect_equal(result, expected)

})

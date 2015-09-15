context("basic HMM functions in AIL")

test_that("AIL check_geno works", {

    # autosome
    for(i in 0:5)
        expect_true(test_check_geno("ail", i, TRUE, FALSE, FALSE, c(20, 0)))
    for(i in 1:3)
        expect_true(test_check_geno("ail", i, FALSE, FALSE, FALSE, c(20, 0)))
    for(i in c(-1, 6))
        expect_false(test_check_geno("ail", i, TRUE, FALSE, FALSE, c(20, 0)))
    for(i in c(0, 4))
        expect_false(test_check_geno("ail", i, FALSE, FALSE, FALSE, c(20, 0)))

    for(dir in 0:2) {
        # X chromosome female
        for(i in 0:5)
            expect_true(test_check_geno("ail", i, TRUE, TRUE, TRUE, c(20, dir)))
        for(i in 1:3)
            expect_true(test_check_geno("ail", i, FALSE, TRUE, TRUE, c(20, dir)))
        for(i in c(-1, 6))
            expect_false(test_check_geno("ail", i, TRUE, TRUE, TRUE, c(20, dir)))
        for(i in c(0, 4:6))
            expect_false(test_check_geno("ail", i, FALSE, TRUE, TRUE, c(20, dir)))

        # X chromosome male
        for(i in 0:5)
            expect_true(test_check_geno("ail", i, TRUE, TRUE, FALSE, c(20, dir)))
        for(i in c(4,5))
            expect_true(test_check_geno("ail", i, FALSE, TRUE, FALSE, c(20, dir)))
        for(i in c(-1, 6))
            expect_false(test_check_geno("ail", i, TRUE, TRUE, FALSE, c(20, dir)))
        for(i in c(0:3, 6))
            expect_false(test_check_geno("ail", i, FALSE, TRUE, FALSE, c(20, dir)))
    }

})

test_that("AIL n_gen works", {

    expect_equal(test_ngen("ail", FALSE), 3)
    expect_equal(test_ngen("ail", TRUE),  5)

})

test_that("AIL possible_gen works", {

    # autosome
    expect_equal(test_possible_gen("ail", FALSE, FALSE, c(20,0)), 1:3)

    for(dir in 0:2) {
        # X female
        expect_equal(test_possible_gen("ail", TRUE, TRUE, c(20,dir)), 1:3)
        # X male
        expect_equal(test_possible_gen("ail", TRUE, FALSE, c(20,dir)), 4:5)
    }

})
test_that("AIL init works", {

    for(n_gen in c(3,9,12,15)) {
        # autosome
        expect_equal(test_init("ail", 1, FALSE, FALSE, c(n_gen, 0)), log(0.25))
        expect_equal(test_init("ail", 2, FALSE, FALSE, c(n_gen, 0)), log(0.5))
        expect_equal(test_init("ail", 3, FALSE, FALSE, c(n_gen, 0)), log(0.25))

        # X chr balanced
        #     female
        expect_equal(test_init("ail", 1, TRUE, TRUE, c(n_gen, 2)), log(0.25))
        expect_equal(test_init("ail", 2, TRUE, TRUE, c(n_gen, 2)), log(0.5))
        expect_equal(test_init("ail", 3, TRUE, TRUE, c(n_gen, 2)), log(0.25))

        #     male
        expect_equal(test_init("ail", 4, TRUE, FALSE, c(n_gen, 2)), log(0.5))
        expect_equal(test_init("ail", 5, TRUE, FALSE, c(n_gen, 2)), log(0.5))
    }

    for(n_gen in c(3,9,12,15)) {
        fprob <- (2/3) + (1/3)*(-1/2)^n_gen
        mprob <- (2/3) + (1/3)*(-1/2)^(n_gen-1)

        # X chr, AxB
        #     female
        expect_equal(test_init("ail", 1, TRUE, TRUE, c(n_gen, 0)), log(fprob^2))
        expect_equal(test_init("ail", 2, TRUE, TRUE, c(n_gen, 0)), log(2*fprob*(1-fprob)))
        expect_equal(test_init("ail", 3, TRUE, TRUE, c(n_gen, 0)), log((1-fprob)^2))

        #     male
        expect_equal(test_init("ail", 4, TRUE, FALSE, c(n_gen, 0)), log(mprob))
        expect_equal(test_init("ail", 5, TRUE, FALSE, c(n_gen, 0)), log(1-mprob))

        # X chr, BxA
        #     female
        expect_equal(test_init("ail", 3, TRUE, TRUE, c(n_gen, 1)), log(fprob^2))
        expect_equal(test_init("ail", 2, TRUE, TRUE, c(n_gen, 1)), log(2*fprob*(1-fprob)))
        expect_equal(test_init("ail", 1, TRUE, TRUE, c(n_gen, 1)), log((1-fprob)^2))

        #     male
        expect_equal(test_init("ail", 5, TRUE, FALSE, c(n_gen, 1)), log(mprob))
        expect_equal(test_init("ail", 4, TRUE, FALSE, c(n_gen, 1)), log(1-mprob))
    }

})


test_that("AIL emit works", {

    # autosome
    eps <- 0.01
    for(i in 1:3)
        expect_equal(test_emit("ail", 0, i, eps, integer(0), FALSE, FALSE, c(20,0)), 0)

    expect_equal(test_emit("ail", 1, 1, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 2, 1, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 3, 1, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 4, 1, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps/2))
    expect_equal(test_emit("ail", 5, 1, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps))

    expect_equal(test_emit("ail", 1, 2, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 2, 2, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 3, 2, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 4, 2, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps/2))
    expect_equal(test_emit("ail", 5, 2, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps/2))

    expect_equal(test_emit("ail", 1, 3, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 2, 3, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 3, 3, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 4, 3, eps, integer(0), FALSE, FALSE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 5, 3, eps, integer(0), FALSE, FALSE, c(20,0)), log(1-eps/2))


    # X female
    for(i in 1:3)
        expect_equal(test_emit("ail", 0, i, eps, integer(0), TRUE, TRUE, c(20,0)), 0)

    expect_equal(test_emit("ail", 1, 1, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 2, 1, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 3, 1, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 4, 1, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps/2))
    expect_equal(test_emit("ail", 5, 1, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps))

    expect_equal(test_emit("ail", 1, 2, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 2, 2, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 3, 2, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 4, 2, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps/2))
    expect_equal(test_emit("ail", 5, 2, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps/2))

    expect_equal(test_emit("ail", 1, 3, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 2, 3, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps/2))
    expect_equal(test_emit("ail", 3, 3, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 4, 3, eps, integer(0), TRUE, TRUE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 5, 3, eps, integer(0), TRUE, TRUE, c(20,0)), log(1-eps/2))

    # X male
    for(i in 4:5)
        expect_equal(test_emit("ail", 0, i, eps, integer(0), TRUE, FALSE, c(20,0)), 0)
    expect_equal(test_emit("ail", 1, 4, eps, integer(0), TRUE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 2, 4, eps, integer(0), TRUE, FALSE, c(20,0)), 0)
    expect_equal(test_emit("ail", 3, 4, eps, integer(0), TRUE, FALSE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 4, 4, eps, integer(0), TRUE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 5, 4, eps, integer(0), TRUE, FALSE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 1, 5, eps, integer(0), TRUE, FALSE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 2, 5, eps, integer(0), TRUE, FALSE, c(20,0)), 0)
    expect_equal(test_emit("ail", 3, 5, eps, integer(0), TRUE, FALSE, c(20,0)), log(1-eps))
    expect_equal(test_emit("ail", 4, 5, eps, integer(0), TRUE, FALSE, c(20,0)), log(eps))
    expect_equal(test_emit("ail", 5, 5, eps, integer(0), TRUE, FALSE, c(20,0)), log(1-eps))

})

test_that("AIL step works", {

    # autosome
    for(rf in c(0.01, 0.0001)) {
        for(ngen in c(3, 9, 12, 15)) {
            R <- 1 - 2*( 0.25*(1 + (1-2*rf)*(1-rf)^(ngen-2)))

            expect_equal(test_step("ail", 1, 1, rf, FALSE, FALSE, c(ngen, 0)), log((1-R)^2))
            expect_equal(test_step("ail", 1, 2, rf, FALSE, FALSE, c(ngen, 0)), log(2*R*(1-R)))
            expect_equal(test_step("ail", 1, 3, rf, FALSE, FALSE, c(ngen, 0)), log(R^2))
            expect_equal(test_step("ail", 2, 1, rf, FALSE, FALSE, c(ngen, 0)), log(R*(1-R)))
            expect_equal(test_step("ail", 2, 2, rf, FALSE, FALSE, c(ngen, 0)), log((1-R)^2+R^2))
            expect_equal(test_step("ail", 2, 3, rf, FALSE, FALSE, c(ngen, 0)), log(R*(1-R)))
            expect_equal(test_step("ail", 3, 1, rf, FALSE, FALSE, c(ngen, 0)), log(R^2))
            expect_equal(test_step("ail", 3, 2, rf, FALSE, FALSE, c(ngen, 0)), log(2*R*(1-R)))
            expect_equal(test_step("ail", 3, 3, rf, FALSE, FALSE, c(ngen, 0)), log((1-R)^2))
        }
    }

    # X chromosome, balanced
    for(rf in c(0.01, 0.0001)) {
        for(ngen in c(3, 9, 12, 15)) {
            z <- sqrt((1-rf)*(9-rf))
            w <- (1 - rf + z)/4
            y <- (1 - rf - z)/4
            mR <- 1 - 0.25*(2 + (1-2*rf)*(w^(ngen-2) + y^(ngen-2)) +
                            (3 - 5*rf + 2*rf^2)/z*(w^(ngen-2) - y^(ngen-2)))
            fR <- 1 - 0.25*(2 + (1-2*rf)*(w^(ngen-2) + y^(ngen-2)) +
                            (3 - 6*rf + rf^2)/z*(w^(ngen-2) - y^(ngen-2)))

            expect_equal(test_step("ail", 1, 1, rf, TRUE, TRUE, c(ngen, 2)), log((1-fR)^2))
            expect_equal(test_step("ail", 1, 2, rf, TRUE, TRUE, c(ngen, 2)), log(2*fR*(1-fR)))
            expect_equal(test_step("ail", 1, 3, rf, TRUE, TRUE, c(ngen, 2)), log(fR^2))
            expect_equal(test_step("ail", 2, 1, rf, TRUE, TRUE, c(ngen, 2)), log(fR*(1-fR)))
            expect_equal(test_step("ail", 2, 2, rf, TRUE, TRUE, c(ngen, 2)), log((1-fR)^2+fR^2))
            expect_equal(test_step("ail", 2, 3, rf, TRUE, TRUE, c(ngen, 2)), log(fR*(1-fR)))
            expect_equal(test_step("ail", 3, 1, rf, TRUE, TRUE, c(ngen, 2)), log(fR^2))
            expect_equal(test_step("ail", 3, 2, rf, TRUE, TRUE, c(ngen, 2)), log(2*fR*(1-fR)))
            expect_equal(test_step("ail", 3, 3, rf, TRUE, TRUE, c(ngen, 2)), log((1-fR)^2))

            expect_equal(test_step("ail", 4, 4, rf, TRUE, FALSE, c(ngen, 2)), log(1-mR))
            expect_equal(test_step("ail", 4, 5, rf, TRUE, FALSE, c(ngen, 2)), log(mR))
            expect_equal(test_step("ail", 5, 4, rf, TRUE, FALSE, c(ngen, 2)), log(mR))
            expect_equal(test_step("ail", 5, 5, rf, TRUE, FALSE, c(ngen, 2)), log(1-mR))

        }
    }


    # calculate Pr(A) for females, AxB
    calc_q <- function(n_gen)
        2/3 + (1/3)*(-1/2)^n_gen

    # calc Pr(AA haplotype) in males and females, recursively
    calc_p11 <- function(n_gen, rf) {
        if(n_gen == 1) return(c(1, 0.5))
        last <- calc_p11(n_gen-1, rf)
        z <- c((1-rf)*last[2] + rf*calc_q(n_gen-2)*calc_q(n_gen-3),
               0.5*last[1] + (1-rf)/2*last[2] + (rf/2)*calc_q(n_gen-2)*calc_q(n_gen-3))
        z
    }

    # X chromosome, AxB and BxA
    for(rf in c(0.01, 0.0001)) {
        for(ngen in c(3, 9, 12, 15)) {

            qf <- calc_q(ngen)
            qm <- calc_q(ngen-1) # Pr(A) in males is Pr(A) in females for prev generation
            p11 <- calc_p11(ngen, rf) # Pr(AA haplotype) in males and females
            m11 <- p11[1]
            f11 <- p11[2]

            # turn into conditional probabilities along one haplotype
            m1to1 <- m11/qm
            m1to2 <- 1 - m1to1
            m2to1 <- (qm - m11)/(1-qm)
            m2to2 <- 1 - m2to1

            f1to1 <- f11/qf
            f1to2 <- 1 - f1to1
            f2to1 <- (qf - f11)/(1-qf)
            f2to2 <- 1 - f2to1

            expect_equal(test_step("ail", 1, 1, rf, TRUE, TRUE, c(ngen, 0)), log(f1to1^2))
            expect_equal(test_step("ail", 1, 2, rf, TRUE, TRUE, c(ngen, 0)), log(2*f1to1*f1to2))
            expect_equal(test_step("ail", 1, 3, rf, TRUE, TRUE, c(ngen, 0)), log(f1to2^2))
            expect_equal(test_step("ail", 2, 1, rf, TRUE, TRUE, c(ngen, 0)), log(f1to1*f2to1))
            expect_equal(test_step("ail", 2, 2, rf, TRUE, TRUE, c(ngen, 0)), log(f2to1*f1to2 + f1to1*f2to2))
            expect_equal(test_step("ail", 2, 3, rf, TRUE, TRUE, c(ngen, 0)), log(f1to2*f2to2))
            expect_equal(test_step("ail", 3, 1, rf, TRUE, TRUE, c(ngen, 0)), log(f2to1^2))
            expect_equal(test_step("ail", 3, 2, rf, TRUE, TRUE, c(ngen, 0)), log(2*f2to2*f2to1))
            expect_equal(test_step("ail", 3, 3, rf, TRUE, TRUE, c(ngen, 0)), log(f2to2^2))

            expect_equal(test_step("ail", 4, 4, rf, TRUE, FALSE, c(ngen, 0)), log(m1to1))
            expect_equal(test_step("ail", 4, 5, rf, TRUE, FALSE, c(ngen, 0)), log(m1to2))
            expect_equal(test_step("ail", 5, 4, rf, TRUE, FALSE, c(ngen, 0)), log(m2to1))
            expect_equal(test_step("ail", 5, 5, rf, TRUE, FALSE, c(ngen, 0)), log(m2to2))

            expect_equal(test_step("ail", 1, 1, rf, TRUE, TRUE, c(ngen, 1)), log(f2to2^2))
            expect_equal(test_step("ail", 1, 2, rf, TRUE, TRUE, c(ngen, 1)), log(2*f2to2*f2to1))
            expect_equal(test_step("ail", 1, 3, rf, TRUE, TRUE, c(ngen, 1)), log(f2to1^2))
            expect_equal(test_step("ail", 2, 1, rf, TRUE, TRUE, c(ngen, 1)), log(f2to2*f1to2))
            expect_equal(test_step("ail", 2, 2, rf, TRUE, TRUE, c(ngen, 1)), log(f1to2*f2to1 + f1to1*f2to2))
            expect_equal(test_step("ail", 2, 3, rf, TRUE, TRUE, c(ngen, 1)), log(f2to1*f1to1))
            expect_equal(test_step("ail", 3, 1, rf, TRUE, TRUE, c(ngen, 1)), log(f1to2^2))
            expect_equal(test_step("ail", 3, 2, rf, TRUE, TRUE, c(ngen, 1)), log(2*f1to1*f1to2))
            expect_equal(test_step("ail", 3, 3, rf, TRUE, TRUE, c(ngen, 1)), log(f1to1^2))

            expect_equal(test_step("ail", 4, 4, rf, TRUE, FALSE, c(ngen, 1)), log(m2to2))
            expect_equal(test_step("ail", 4, 5, rf, TRUE, FALSE, c(ngen, 1)), log(m2to1))
            expect_equal(test_step("ail", 5, 4, rf, TRUE, FALSE, c(ngen, 1)), log(m1to2))
            expect_equal(test_step("ail", 5, 5, rf, TRUE, FALSE, c(ngen, 1)), log(m1to1))

        }
    }


})

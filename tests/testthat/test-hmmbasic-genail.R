context("basic HMM functions in general AIL")

### TODO none of this has been revised yet

test_that("genail n_gen, n_alleles work", {

    expect_equal(nalleles("genail38"), 38)
    expect_equal(nalleles("genail6"), 6)

    expect_equal(test_ngen("genail38", FALSE), 741)
    expect_equal(test_ngen("genail8", FALSE), 36)
    expect_equal(test_ngen("genail38", TRUE), 779)
    expect_equal(test_ngen("genail8", TRUE), 44)

    # throw error for "genail1"
    expect_error(nalleles("genail1"))

})

test_that("genail possible_gen work", {

    expect_equal(test_possible_gen("genail38", FALSE, FALSE, c(15,rep(1,38))), 1:741)
    expect_equal(test_possible_gen("genail38", TRUE, TRUE, c(15,rep(1,38))), 1:741)
    expect_equal(test_possible_gen("genail38", TRUE, FALSE, c(15,rep(1,38))), 742:779)

})

# FIX_ME: test check_geno

test_that("genail init work", {

    set.seed(20181105)
    alpha <- sample(1:10, 38, replace=TRUE)
    init<- alpha/sum(alpha)

    calcFX <- calcA <- expected <- matrix(nrow=38, ncol=38)
    calcMX <- rep(NA, 38)

    for(i in 1:38) {
        for(j in i:38) {
            g <- mpp_encode_alleles(i, j, 38, FALSE)

            calcA[i,j] <- test_init("genail38", g, FALSE, FALSE, c(8, alpha))
            expected[i,j] <- ifelse(i==j, 2*log(init[i]), log(2) + log(init[i]) + log(init[j]))

            calcFX[i,j] <- test_init("genail38", g, TRUE, TRUE, c(8, alpha))

        }

        calcMX[i] <- test_init("genail38", i+741, TRUE, FALSE, c(8, alpha))
    }

    expect_equal(calcMX, log(init))
    expect_equal(calcA, expected)
    expect_equal(calcFX, expected)

})

# FIX_ME: test emit

test_that("genail step works", {

    nf <- 7
    ng <- nf + choose(nf,2)
    alpha_int <- sample(1:10, nf, replace=TRUE)
    alpha <- alpha_int/sum(alpha_int)

    g <- t(sapply(1:ng, mpp_decode_geno, nf, FALSE)) # genotypes

    for(rf in c(0.01, 0.1, 0.45)) {
        for(ngen in c(3, 5)) {

            # male X chr
            R <- alpha*(1-(1-rf)^ngen)
            expected <- matrix(rep(alpha*(1-(1-rf)^ngen),nf), ncol=nf, byrow=TRUE)
            diag(expected) <- alpha + (1-alpha)*(1-rf)^ngen

            result <- matrix(ncol=nf, nrow=nf)
            for(i in 1:nf) {
                for(j in 1:nf) {
                    result[i,j] <- test_step(paste0("genail", nf), i+ng, j+ng, rf, TRUE, FALSE, c(ngen, alpha_int))
                }
            }
            # rows sum to 1?
            expect_equal(rowSums(exp(result)), rep(1, nf))

            # matches what I expected?
            expect_equal(result, log(expected))


            # autosome
            result <- matrix(ncol=ng, nrow=ng)
            for(i in 1:ng) {
                for(j in 1:ng) {
                    result[i,j] <- test_step(paste0("genail", nf), i, j, rf, FALSE, FALSE, c(ngen, alpha_int))
                }
            }
            # rows sum to 1?
            expect_equal(rowSums(exp(result)), rep(1, ng))
        }
    }

})


test_that("genail geno_names work", {

    # if 38 founders, using upper case and lower case letters. Ugh.
    alleles <- c(LETTERS, letters[1:(38-26)])
    expected <- sapply(1:741, function(a) paste(alleles[mpp_decode_geno(a, 38, FALSE)], collapse=""))
    expect_equal( geno_names("genail38", alleles, FALSE), expected)

    hemi <- paste(alleles, "Y", sep="")
    expect_equal( geno_names("genail38", alleles, TRUE), c(expected, hemi))

    # could also use two-letter allele codes, but ugly
    alleles <- c(paste0("A",LETTERS), paste0("B", LETTERS[1:(38-26)]))
    expected <- sapply(1:741, function(a) paste(alleles[mpp_decode_geno(a, 38, FALSE)], collapse=""))
    expect_equal( geno_names("genail38", alleles, FALSE), expected)

    hemi <- paste(alleles, "Y", sep="")
    expect_equal( geno_names("genail38", alleles, TRUE), c(expected, hemi))

})

test_that("genail nrec work", {

    nf <- 19
    ng <- nf + choose(nf, 2)
    crosstype <- paste0("genail", nf)

    # male X chromosome
    calculated <- matrix(nrow=nf, ncol=nf)
    for(i in 1:nf) {
        for(j in 1:nf) {
            calculated[i,j] <- test_nrec(crosstype, i+ng, j+ng, TRUE, FALSE, c(20, rep(1, nf)))
        }
    }
    expect_equal(calculated, 1-diag(nf))

    calculatedA <- calculatedfX <- expected <- matrix(nrow=ng, ncol=ng)
    g <- sapply(1:ng, function(g) mpp_decode_geno(g, nf, FALSE))
    # female X chromosome or autosome
    for(i in 1:ng) {
        for(j in 1:ng) {

            calculatedA[i,j] <- test_nrec(crosstype, i, j, FALSE, FALSE, c(20, rep(1, nf)))
            calculatedfX[i,j] <- test_nrec(crosstype, i, j, TRUE, TRUE, c(20, rep(1, nf)))

            expected[i,j] <- min(sum(g[,i]!=g[,j]), sum(g[,i]!=rev(g[,j])))
        }
    }
    expect_equal(calculatedA, expected)
    expect_equal(calculatedfX, expected)

})

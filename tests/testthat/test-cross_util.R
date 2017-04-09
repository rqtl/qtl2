# test multi-parent population utilities
context("test MPP utilities")

g <- matrix(ncol=8, nrow=8)
cur <- 1
for(i in 1:8) {
    for(j in 1:i) {
        g[j,i] <- g[i,j] <- cur
        cur <- cur + 1
    }
}

test_that("mpp_encode_alleles works for 8 alleles, phase unknown", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_encode_alleles(i, i, 8, FALSE), homg[i])

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_encode_alleles(i, j, 8, FALSE), g[i,j])
        }
    }
})

test_that("mpp_decode_geno works for 8 alleles, phase unknown", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_decode_geno(homg[i], 8, FALSE), c(i,i))

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_decode_geno(g[i, j], 8, FALSE), sort(c(i,j)))
        }
    }
})

test_that("mpp_is_het works for 8 alleles, phase unknown", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_is_het(homg[i], 8, FALSE), FALSE)

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_is_het(g[i, j], 8, FALSE), ifelse(i==j, FALSE, TRUE))
        }
    }
})

test_that("mpp_geno_names works for 8 alleles", {
    a <- LETTERS[1:8]

    # genotype names
    gn <- paste0(a[row(g)], a[col(g)])[upper.tri(g, TRUE)]

    expect_equal(mpp_geno_names(a, FALSE), gn)
    expect_equal(mpp_geno_names(a, TRUE), c(gn, paste0(a, "Y")))
})


# phase-known case
cur <- 37
for(i in 2:8) {
    for(j in 1:(i-1)) {
        g[i,j] <- cur
        cur <- cur + 1
    }
}

test_that("mpp_encode_alleles works for 8 alleles, phase known", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_encode_alleles(i, i, 8, TRUE), homg[i])

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_encode_alleles(i, j, 8, TRUE), g[i,j])
        }
    }
})

test_that("mpp_decode_geno works for 8 alleles, phase known", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_decode_geno(homg[i], 8, TRUE), c(i,i))

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_decode_geno(g[i, j], 8, TRUE), c(i,j))
        }
    }
})

test_that("mpp_is_het works for 8 alleles, phase known", {
    homg <- c(1,3,6,10,15,21,28,36)
    for(i in seq(along=homg))
        expect_equal(mpp_is_het(homg[i], 8, TRUE), FALSE)

    for(i in 1:8) {
        for(j in 1:8) {
            expect_equal(mpp_is_het(g[i, j], 8, TRUE), ifelse(i==j, FALSE, TRUE))
        }
    }
})

test_that("invert_founder_index works", {

    expect_equal(invert_founder_index(c(2,3,1,4)), c(2,0,1,3))
    expect_equal(invert_founder_index(c(2,4,3,1)), c(3,0,2,1))
    expect_equal(invert_founder_index(c(7,8,3,5,4,1,6,2)), c(5,7,2,4,3,6,0,1))

})

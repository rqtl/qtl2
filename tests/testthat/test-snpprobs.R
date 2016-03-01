context("genoprobs -> snpprobs")

test_that("calc_sdp works", {

    expect_error(calc_sdp(rbind(0,1,0))) # should be at least two columns

    expect_equal(calc_sdp( rbind( c(0,0), c(1,0), c(0,1), c(1,1) ) ),
                 c(0, 1, 2, 3) )

    expect_equal(calc_sdp( rbind(c(1, 0, 0, 0),
                                 c(0, 1, 0, 0),
                                 c(0, 0, 1, 0),
                                 c(0, 0, 0, 1),
                                 c(1, 1, 0, 0),
                                 c(1, 0, 1, 0),
                                 c(1, 0, 0, 1),
                                 c(0, 1, 1, 0),
                                 c(0, 1, 0, 1),
                                 c(0, 0, 1, 1),
                                 c(1, 1, 1, 0),
                                 c(1, 1, 0, 1),
                                 c(1, 0, 1, 1),
                                 c(0, 1, 1, 1),
                                 c(1, 1, 1, 1)) ),
                 c(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14, 15) )

    g <- rbind(c(1,0,0,0,0,0,0,0), # 1
               c(0,1,0,0,0,0,0,0), # 2
               c(0,0,0,0,0,1,0,0), # 32
               c(0,0,0,0,0,0,0,1), # 128
               c(1,0,0,0,0,0,0,1), # 129
               c(0,1,0,1,0,1,0,1), # 170
               c(1,0,1,0,1,0,1,0), # 85
               c(1,1,1,1,0,0,0,0), # 15
               c(0,0,0,0,1,1,1,1), # 240
               c(1,1,0,0,0,0,1,1)) # 195

    expect_equal(calc_sdp(g), c(1,2,32,128,129,170,85,15,240,195) )

    set.seed(38444584)
    g <- matrix(sample(0:1, 8*12, replace=TRUE), ncol=8)
    expect_equal(calc_sdp(g),
                 apply(g, 1, function(a) sum(a*2^(seq(along=a)-1))))
})


test_that(".alleleprob_to_snpprob works in simple cases", {

    # conversion correct, with 1 map position
    pr <- array(c(0.15, 0.03, 0.20, 0.08, 0.18, 0.09, 0.15, 0.11),
                dim=c(1, 8, 1))
    sdp <- 1:255
    interval <- rep(0, length(sdp))
    on_map <- rep(TRUE, length(sdp))

    # sdp -> genotype matrix
    g <- matrix(rep(sdp,rep(8, length(sdp))), ncol=length(sdp))
    g <- (g %/% 2^(0:7)) %% 2

    expected <- array(dim=c(1, 2, length(sdp)))
    for(i in 1:length(sdp)) {
        expected[1,1,i] <- sum(pr[1,g[,i]==0,1])
        expected[1,2,i] <- sum(pr[1,g[,i]==1,1])
    }

    expect_equal(.alleleprob_to_snpprob(pr, sdp, interval, on_map),
                 expected)


    # add a second map position
    pr <- array(c(0.15, 0.03, 0.20, 0.08, 0.18, 0.09, 0.15, 0.11,
                  0.07, 0.05, 0.17, 0.02, 0.46, 0.03, 0.11, 0.08),
                dim=c(1, 8, 2))

    prsum <- (pr[,,1] + pr[,,2])/2
    interval[201:255] <- 1   # last 50 at 2nd marker
    on_map[151:200] <- FALSE # 50 before that in-between

    expected <- array(dim=c(1, 2, length(sdp)))
    for(i in 1:150) {
        expected[1,1,i] <- sum(pr[1,g[,i]==0,1])
        expected[1,2,i] <- sum(pr[1,g[,i]==1,1])
    }
    for(i in 151:200) {
        expected[1,1,i] <- sum((pr[1,g[,i]==0,1] + pr[1,g[,i]==0,2])/2)
        expected[1,2,i] <- sum((pr[1,g[,i]==1,1] + pr[1,g[,i]==1,2])/2)
    }
    for(i in 201:255) {
        expected[1,1,i] <- sum(pr[1,g[,i]==0,2])
        expected[1,2,i] <- sum(pr[1,g[,i]==1,2])
    }

    expect_equal(.alleleprob_to_snpprob(pr, sdp, interval, on_map),
                 expected)

})

test_that("genocol_to_snpcol works with 8 alleles", {

    # A has one allele and rest have another
    expected <- rep(0, 36)
    expected[1] <- 2
    expected[c(2,4,7,11,16,22,29)] <- 1
    expect_equal(genocol_to_snpcol(8,1), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,254), 2-expected)

    # B has one allele and rest have another
    expected <- rep(0, 36)
    expected[3] <- 2
    expected[c(2,5,8,12,17,23,30)] <- 1
    expect_equal(genocol_to_snpcol(8,2), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,253), 2-expected)

    # E has one allele and rest have another
    expected <- rep(0, 36)
    expected[15] <- 2
    expected[c(11,12,13,14,20,26,33)] <- 1
    expect_equal(genocol_to_snpcol(8,16), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,239), 2-expected)

    # H has one allele and rest have another
    expected <- rep(0, 36)
    expected[36] <- 2
    expected[29:35] <- 1
    expect_equal(genocol_to_snpcol(8,128), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,127), 2-expected)

    # ABCD have one allele and EFGH has
    expected <- rep(1, 36)
    expected[1:10] <- 2
    expected[c(15,20,21,26:28,33:36)] <- 0
    expect_equal(genocol_to_snpcol(8,15), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,255-15), 2-expected)

    # ACEG have one allele and BDFH has
    expected <- rep(1, 36)
    expected[c(1,4,6,11,13,15,22,24,26,28)] <- 2
    expected[c(3,8,10,17,19,21,30,32,34,36)] <- 0
    expect_equal(genocol_to_snpcol(8,85), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,255-85), 2-expected)

    # BE have one allele and ACDFGH have other
    expected <- rep(0, 36)
    expected[c(3,12,15)] <- 2
    expected[c(2,5,8,17,23,30,11,13,14,20,26,33)] <- 1
    expect_equal(genocol_to_snpcol(8,18), expected)
    # opposite
    expect_equal(genocol_to_snpcol(8,255-18), 2-expected)

    # test all opposites
    for(g in 1:254)
        expect_equal(genocol_to_snpcol(8, g),
                     2 - genocol_to_snpcol(8, 255-g))

})

test_that("genocol_to_snpcol works with 4 alleles", {

    # A,BCD (1)
    expect_equal(genocol_to_snpcol(4, 1),
                 c(2,1,0,1,0,0,1,0,0,0))
    expect_equal(genocol_to_snpcol(4, 15-1),
                 2-c(2,1,0,1,0,0,1,0,0,0))
    # B,ACD (2)
    expect_equal(genocol_to_snpcol(4, 2),
                 c(0,1,2,0,1,0,0,1,0,0))
    expect_equal(genocol_to_snpcol(4, 15-2),
                 2-c(0,1,2,0,1,0,0,1,0,0))
    # C,ABD (4)
    expect_equal(genocol_to_snpcol(4, 4),
                 c(0,0,0,1,1,2,0,0,1,0))
    expect_equal(genocol_to_snpcol(4, 15-4),
                 2-c(0,0,0,1,1,2,0,0,1,0))
    # D,ABC (8)
    expect_equal(genocol_to_snpcol(4, 8),
                 c(0,0,0,0,0,0,1,1,1,2))
    expect_equal(genocol_to_snpcol(4, 15-8),
                 2-c(0,0,0,0,0,0,1,1,1,2))
    # AB,CD (3)
    expect_equal(genocol_to_snpcol(4, 3),
                 c(2,2,2,1,1,0,1,1,0,0))
    expect_equal(genocol_to_snpcol(4, 15-3),
                 2-c(2,2,2,1,1,0,1,1,0,0))
    # AC,BD (5)
    expect_equal(genocol_to_snpcol(4, 5),
                 c(2,1,0,2,1,2,1,0,1,0))
    expect_equal(genocol_to_snpcol(4, 15-5),
                 2-c(2,1,0,2,1,2,1,0,1,0))
    # AD,BC (9)
    expect_equal(genocol_to_snpcol(4, 9),
                 c(2,1,0,1,0,0,2,1,1,2))
    expect_equal(genocol_to_snpcol(4, 15-9),
                 2-c(2,1,0,1,0,0,2,1,1,2))
})

test_that(".genoprob_to_snpprob works with 8 alleles", {

    # make a fake probability matrix
    set.seed(36319176)
    npos <- 8
    nind <- 5
    x <- matrix(runif(npos*nind*36), ncol=36)
    x <- x/rowSums(x)
    pr <- array(dim=c(nind, 36, npos))
    k <- 1
    for(i in 1:nind) {
        for(j in 1:npos) {
            pr[i,,j] <- x[k,]
            k <- k+1
        }
    }

    # errors
    expect_error(.genoprob_to_snpprob(pr[,-1,], 1, 0, TRUE)) # not enough columns
    expect_error(.genoprob_to_snpprob(pr[,-(1:2),], 1, 0, TRUE)) # not enough columns
    expect_error(.genoprob_to_snpprob(pr, c(1,0,200), c(0,1,2), c(TRUE,FALSE,FALSE))) # SDP out of range
    expect_error(.genoprob_to_snpprob(pr, c(1,2,256), c(0,1,2), c(TRUE,FALSE,FALSE))) # SDP out of range
    expect_error(.genoprob_to_snpprob(pr, c(1,255,200), c(0,1,-1), c(TRUE,FALSE,FALSE))) # snp out of map range
    expect_error(.genoprob_to_snpprob(pr, c(1,2,254), c(0,1,npos-1), c(TRUE,FALSE,FALSE)))    # snp out of map range

    nsnp <- 200
    sdp <- sample(255, nsnp, replace=TRUE)
    interval <- sample(0:(npos-1), nsnp, replace=TRUE)
    on_map <- sample(c(FALSE, TRUE), nsnp, replace=TRUE)
    on_map[interval==npos-1] <- TRUE

    # calculate snp genotype probabilities
    result <- .genoprob_to_snpprob(pr, sdp, interval, on_map)
    expect_equal(dim(result), c(nind,3,nsnp))

    # calculate in R
    expected <- array(dim=c(nind, 3, nsnp))
    for(snp in 1:nsnp) {
        snpcol <- genocol_to_snpcol(8, sdp[snp])
        for(snpg in 1:3) {
            if(on_map[snp])
                expected[,snpg,snp] <- rowSums(pr[,snpcol==snpg-1,interval[snp]+1,drop=FALSE])
            else
                expected[,snpg,snp] <- rowSums((pr[,snpcol==snpg-1,interval[snp]+1,drop=FALSE] +
                                                pr[,snpcol==snpg-1,interval[snp]+2,drop=FALSE])/2)
        }
    }
    expect_equal(result, expected)

})


test_that("Xgenocol_to_snpcol works with 8 alleles", {

    # flip expected when flipping SDP
    opp <- function(x) c(2-x[1:36],7-x[37:44])

    # A has one allele and rest have another
    expected <- rep(0, 44)
    expected[1] <- 2
    expected[c(2,4,7,11,16,22,29)] <- 1
    expected[37:44] <- c(4,3,3,3,3,3,3,3)
    expect_equal(Xgenocol_to_snpcol(8,1), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,254), opp(expected))

    # B has one allele and rest have another
    expected <- rep(0, 44)
    expected[3] <- 2
    expected[c(2,5,8,12,17,23,30)] <- 1
    expected[37:44] <- c(3,4,3,3,3,3,3,3)
    expect_equal(Xgenocol_to_snpcol(8,2), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,253), opp(expected))

    # E has one allele and rest have another
    expected <- rep(0, 44)
    expected[15] <- 2
    expected[c(11,12,13,14,20,26,33)] <- 1
    expected[37:44] <- c(3,3,3,3,4,3,3,3)
    expect_equal(Xgenocol_to_snpcol(8,16), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,239), opp(expected))

    # H has one allele and rest have another
    expected <- rep(0, 44)
    expected[36] <- 2
    expected[29:35] <- 1
    expected[37:44] <- c(3,3,3,3,3,3,3,4)
    expect_equal(Xgenocol_to_snpcol(8,128), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,127), opp(expected))

    # ABCD have one allele and EFGH has
    expected <- rep(1, 44)
    expected[1:10] <- 2
    expected[c(15,20,21,26:28,33:36)] <- 0
    expected[37:44] <- c(4,4,4,4,3,3,3,3)
    expect_equal(Xgenocol_to_snpcol(8,15), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,255-15), opp(expected))

    # ACEG have one allele and BDFH has
    expected <- rep(1, 44)
    expected[c(1,4,6,11,13,15,22,24,26,28)] <- 2
    expected[c(3,8,10,17,19,21,30,32,34,36)] <- 0
    expected[37:44] <- c(4,3,4,3,4,3,4,3)
    expect_equal(Xgenocol_to_snpcol(8,85), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,255-85), opp(expected))

    # BE have one allele and ACDFGH have other
    expected <- rep(0, 36)
    expected[c(3,12,15)] <- 2
    expected[c(2,5,8,17,23,30,11,13,14,20,26,33)] <- 1
    expected[37:44] <- c(3,4,3,3,4,3,3,3)
    expect_equal(Xgenocol_to_snpcol(8,18), expected)
    # opposite
    expect_equal(Xgenocol_to_snpcol(8,255-18), opp(expected))

    # test all opposites
    for(g in 1:254)
        expect_equal(Xgenocol_to_snpcol(8, g),
                     opp(Xgenocol_to_snpcol(8, 255-g)))

})

test_that("Xgenocol_to_snpcol works with 4 alleles", {

    # A,BCD (1)
    expect_equal(Xgenocol_to_snpcol(4, 1),
                 c(2,1,0,1,0,0,1,0,0,0,4,3,3,3))
    expect_equal(Xgenocol_to_snpcol(4, 15-1),
                 c(0,1,2,1,2,2,1,2,2,2,3,4,4,4))
    # B,ACD (2)
    expect_equal(Xgenocol_to_snpcol(4, 2),
                 c(0,1,2,0,1,0,0,1,0,0,3,4,3,3))
    expect_equal(Xgenocol_to_snpcol(4, 15-2),
                 c(2,1,0,2,1,2,2,1,2,2,4,3,4,4))
    # C,ABD (4)
    expect_equal(Xgenocol_to_snpcol(4, 4),
                 c(0,0,0,1,1,2,0,0,1,0,3,3,4,3))
    expect_equal(Xgenocol_to_snpcol(4, 15-4),
                 c(2,2,2,1,1,0,2,2,1,2,4,4,3,4))
    # D,ABC (8)
    expect_equal(Xgenocol_to_snpcol(4, 8),
                 c(0,0,0,0,0,0,1,1,1,2,3,3,3,4))
    expect_equal(Xgenocol_to_snpcol(4, 15-8),
                 c(2,2,2,2,2,2,1,1,1,0,4,4,4,3))
    # AB,CD (3)
    expect_equal(Xgenocol_to_snpcol(4, 3),
                 c(2,2,2,1,1,0,1,1,0,0,4,4,3,3))
    expect_equal(Xgenocol_to_snpcol(4, 15-3),
                 c(0,0,0,1,1,2,1,1,2,2,3,3,4,4))
    # AC,BD (5)
    expect_equal(Xgenocol_to_snpcol(4, 5),
                 c(2,1,0,2,1,2,1,0,1,0,4,3,4,3))
    expect_equal(Xgenocol_to_snpcol(4, 15-5),
                 c(0,1,2,0,1,0,1,2,1,2,3,4,3,4))
    # AD,BC (9)
    expect_equal(Xgenocol_to_snpcol(4, 9),
                 c(2,1,0,1,0,0,2,1,1,2,4,3,3,4))
    expect_equal(Xgenocol_to_snpcol(4, 15-9),
                 c(0,1,2,1,2,2,0,1,1,0,3,4,4,3))
})


test_that(".Xgenoprob_to_snpprob works with 8 alleles", {

    # make a fake probability matrix
    set.seed(20939472)
    npos <- 8
    nind <- 5
    x <- matrix(runif(npos*nind*36), ncol=36)
    x <- x/rowSums(x)
    x2 <- matrix(runif(npos*nind*8), ncol=8)
    x2 <- x2/rowSums(x2)
    x <- rbind(cbind(x, matrix(0, ncol=8, nrow=nrow(x))),
               cbind(matrix(0, ncol=36, nrow=nrow(x2)), x2))
    nind <- nind*2

    pr <- array(dim=c(nind, 44, npos))
    k <- 1
    for(i in 1:nind) {
        for(j in 1:npos) {
            pr[i,,j] <- x[k,]
            k <- k+1
        }
    }

    # errors
    expect_error(.Xgenoprob_to_snpprob(pr[,-1,], 1, 0, TRUE)) # not enough columns
    expect_error(.Xgenoprob_to_snpprob(pr[,-(1:2),], 1, 0, TRUE)) # not enough columns
    expect_error(.Xgenoprob_to_snpprob(pr, c(1,0,200), c(0,1,2), c(TRUE,FALSE,FALSE))) # SDP out of range
    expect_error(.Xgenoprob_to_snpprob(pr, c(1,2,256), c(0,1,2), c(TRUE,FALSE,FALSE))) # SDP out of range
    expect_error(.Xgenoprob_to_snpprob(pr, c(1,255,200), c(0,1,-1), c(TRUE,FALSE,FALSE))) # snp out of map range
    expect_error(.Xgenoprob_to_snpprob(pr, c(1,2,254), c(0,1,npos-1), c(TRUE,FALSE,FALSE)))    # snp out of map range

    nsnp <- 200
    sdp <- sample(255, nsnp, replace=TRUE)
    interval <- sample(0:(npos-1), nsnp, replace=TRUE)
    on_map <- sample(c(FALSE, TRUE), nsnp, replace=TRUE)
    on_map[interval==npos-1] <- TRUE

    # calculate snp genotype probabilities
    result <- .Xgenoprob_to_snpprob(pr, sdp, interval, on_map)
    expect_equal(dim(result), c(nind,5,nsnp))

    # calculate in R
    expected <- array(dim=c(nind, 5, nsnp))
    for(snp in 1:nsnp) {
        snpcol <- Xgenocol_to_snpcol(8, sdp[snp])
        for(snpg in 1:5) {
            if(on_map[snp])
                expected[,snpg,snp] <- rowSums(pr[,snpcol==snpg-1,interval[snp]+1,drop=FALSE])
            else
                expected[,snpg,snp] <- rowSums((pr[,snpcol==snpg-1,interval[snp]+1,drop=FALSE] +
                                                pr[,snpcol==snpg-1,interval[snp]+2,drop=FALSE])/2)
        }
    }
    expect_equal(result, expected)

})

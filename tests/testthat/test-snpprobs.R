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


test_that("alleleprob_to_snpprob works in simple cases", {

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

    expect_equal(alleleprob_to_snpprob(pr, sdp, interval, on_map),
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

    expect_equal(alleleprob_to_snpprob(pr, sdp, interval, on_map),
                 expected)

})

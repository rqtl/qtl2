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

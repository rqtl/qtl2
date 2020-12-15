context("SDP codes")

test_that("calc_sdp works", {

    expect_error(suppressWarnings(calc_sdp(rbind(1,3,1)))) # should be at least two columns

    expect_equal(suppressWarnings( calc_sdp( rbind( c(1,1), c(3,1), c(1,3), c(3,3) ) ) ),
                 c(1, 2) )

    expect_equal(calc_sdp( rbind(c(3, 1, 1, 1),
                                 c(1, 3, 1, 1),
                                 c(1, 1, 3, 1),
                                 c(1, 1, 1, 3),
                                 c(3, 3, 1, 1),
                                 c(3, 1, 3, 1),
                                 c(3, 1, 1, 3),
                                 c(1, 3, 3, 1),
                                 c(1, 3, 1, 3),
                                 c(1, 1, 3, 3),
                                 c(3, 3, 3, 1),
                                 c(3, 3, 1, 3),
                                 c(3, 1, 3, 3),
                                 c(1, 3, 3, 3)) ),
                 c(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14) )

    g <- rbind(c(3,1,1,1,1,1,1,1), # 1
               c(1,3,1,1,1,1,1,1), # 2
               c(1,1,1,1,1,3,1,1), # 32
               c(1,1,1,1,1,1,1,3), # 128
               c(3,1,1,1,1,1,1,3), # 129
               c(1,3,1,3,1,3,1,3), # 170
               c(3,1,3,1,3,1,3,1), # 85
               c(3,3,3,3,1,1,1,1), # 15
               c(1,1,1,1,3,3,3,3), # 240
               c(3,3,1,1,1,1,3,3)) # 195

    expect_equal(calc_sdp(g), c(1,2,32,128,129,170,85,15,240,195) )

    set.seed(38444584)
    g <- matrix(sample(c(1,3), 8*12, replace=TRUE), ncol=8)
    n_AA <- rowSums(g==1)
    g <- g[n_AA > 0 & n_AA < 8,]
    expect_equal(calc_sdp(g),
                 apply(g, 1, function(a) sum(((a-1)/2)*2^(seq(along=a)-1))))

    expect_equal( calc_sdp( c(1,1,1,3,1,1,1,1) ), 8)
    expect_equal( calc_sdp( data.frame(1,1,1,3,1,1,1,1)),  8)

})


test_that("invert_sdp works", {

    expected <- rbind(c(3,1), c(1,3))
    expect_equal(invert_sdp(c(1,2), 2), expected)

    g <- rbind(c(3, 1, 1, 1),
               c(1, 3, 1, 1),
               c(1, 1, 3, 1),
               c(1, 1, 1, 3),
               c(3, 3, 1, 1),
               c(3, 1, 3, 1),
               c(3, 1, 1, 3),
               c(1, 3, 3, 1),
               c(1, 3, 1, 3),
               c(1, 1, 3, 3),
               c(3, 3, 3, 1),
               c(3, 3, 1, 3),
               c(3, 1, 3, 3),
               c(1, 3, 3, 3))
    expect_equal(invert_sdp(c(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14), 4),
                 g)

    g <- rbind(c(3,1,1,1,1,1,1,1), # 1
               c(1,3,1,1,1,1,1,1), # 2
               c(1,1,1,1,1,3,1,1), # 32
               c(1,1,1,1,1,1,1,3), # 128
               c(3,1,1,1,1,1,1,3), # 129
               c(1,3,1,3,1,3,1,3), # 170
               c(3,1,3,1,3,1,3,1), # 85
               c(3,3,3,3,1,1,1,1), # 15
               c(1,1,1,1,3,3,3,3), # 240
               c(3,3,1,1,1,1,3,3)) # 195
    expect_equal(invert_sdp(c(1,2,32,128,129,170,85,15,240,195), 8),
                 g)

    set.seed(38444584)
    g <- matrix(sample(c(1,3), 8*12, replace=TRUE), ncol=8)
    n_AA <- rowSums(g==1)
    g <- g[n_AA > 0 & n_AA < 8,]
    expect_equal(invert_sdp(calc_sdp(g), 8),
                 g)
})


test_that("sdp2char works", {

    expect_equal(sdp2char(c(1,2), 2), c("A|B", "B|A"))

    expect_equal(sdp2char(c(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14), 4),
                 c("A|BCD", "B|ACD", "C|ABD", "D|ABC",
                   "AB|CD", "AC|BD", "AD|BC", "BC|AD", "BD|AC", "CD|AB",
                   "ABC|D", "ABD|C", "ACD|B", "BCD|A"))

    expect_equal(sdp2char(c(1, 2, 4, 8, 3, 5, 9, 6, 10, 12, 7, 11, 13, 14),
                          strains=LETTERS[1:4]),
                 c("A|BCD", "B|ACD", "C|ABD", "D|ABC",
                   "AB|CD", "AC|BD", "AD|BC", "BC|AD", "BD|AC", "CD|AB",
                   "ABC|D", "ABD|C", "ACD|B", "BCD|A"))

    expect_equal(sdp2char(c(1,2,32,128,129,170,85,15,240,195), 8),
                 c("A|BCDEFGH", "B|ACDEFGH", "F|ABCDEGH", "H|ABCDEFG",
                   "AH|BCDEFG", "BDFH|ACEG", "ACEG|BDFH", "ABCD|EFGH",
                   "EFGH|ABCD", "ABGH|CDEF"))

    expect_equal(sdp2char(c(1,2,32,128,129,170,85,15,240,195),
                          strains=LETTERS[1:8]),
                 c("A|BCDEFGH", "B|ACDEFGH", "F|ABCDEGH", "H|ABCDEFG",
                   "AH|BCDEFG", "BDFH|ACEG", "ACEG|BDFH", "ABCD|EFGH",
                   "EFGH|ABCD", "ABGH|CDEF"))

})

context("Chi-square test for column pairs")

test_that("chisq_colpairs works", {

    set.seed(44552625)
    p <- 5
    z <- matrix(sample(1:2, p*100, replace=TRUE), ncol=p)
    colnames(z) <- LETTERS[1:p]
    z[,2] <- 0

    result <- chisq_colpairs(z)

    expected <- matrix(ncol=p,nrow=p)
    cols <- seq(p)[-2]
    for(i in 1:(p-2)) {
        for(j in (i+1):(p-1)) {
            ii <- cols[i]
            jj <- cols[j]
            suppressWarnings(expected[jj,ii] <- expected[ii,jj] <-
                                 chisq.test(table(z[,ii], z[,jj]), correct=FALSE)$stat)
        }
    }
    dimnames(expected) <- list(LETTERS[1:p], LETTERS[1:p])

    expect_equal(result, expected)

    # throws error if you give a single-column matrix
    expect_error( chisq_colpairs(z[,1,drop=FALSE]) )

    # throws error if you give a vector
    expect_error( chisq_colpairs(z[,1]) )

})

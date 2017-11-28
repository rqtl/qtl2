
context("map functions")

test_that("map functions are giving the invertible", {

    mapfunc <- c("haldane", "kosambi", "c-f", "morgan")

    d <- seq(0, 200, by=10)
    for(m in mapfunc) {
        if(m == "morgan") d <- d[d <= 50]

        expect_equal(imf(mf(d, m), m), d)
    }

    r <- seq(0, 0.5, by=0.5)
    for(m in mapfunc)
        expect_equal(imf(mf(r, m), m), r)

})

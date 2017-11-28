context("x-axis positions for scan1 plot")

test_that("xpos_scan1 works", {

    # single chromosome
    expect_equal(xpos_scan1(list("1"=1:5), thechr="1", thepos=1:5), 1:5)

    m <- sort(runif(100, 0, 100))
    expect_equal(xpos_scan1(list("1"=m), thechr="1", thepos=m), m)

    # two chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85))
    expect_equal(xpos_scan1(map, thechr=c("1", "2"), thepos=c(20, 25)),
                 c(22.5, 97.5))
    expect_equal(xpos_scan1(map, thechr="1", thepos=c(15, 40)),
                 c(17.5, 42.5))
    expect_equal(xpos_scan1(map, thechr=c("1", "2"), thepos=10),
                 c(12.5, 82.5))


    expect_error(xpos_scan1(map, thechr=c("1", "2", "2"), thepos=c(5, 10)))
    expect_error(xpos_scan1(map, thechr=c("1", "2"), thepos=c(5, 10, 15)))
    expect_error(xpos_scan1(map, thechr=c("1", "2", "3"), thepos=c(5, 10, 15)))

})

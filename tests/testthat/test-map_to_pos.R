context("Converting map to x-axis positions")

test_that("map_to_xpos works", {

    # single chromosome
    expect_equal(map_to_xpos(list("1"=1:5)), 0:4)

    m <- sort(runif(100, 0, 100))
    expect_equal(map_to_xpos(list("1"=m)), m-min(m))

    # two chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85))
    expect_equal(map_to_xpos(map, 0),
                 c(0, 10, 40, 40, 60, 120))
    expect_equal(map_to_xpos(map, 5),
                 c(2.5, 12.5, 42.5, 47.5, 67.5, 127.5))

    # three chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85),
                "3"=c(90, 100, 150, 175))
    expect_equal(map_to_xpos(map, 0),
                 c(0, 10, 40, 40, 60, 120, 120, 130, 180, 205))
    expect_equal(map_to_xpos(map, 5),
                 c(2.5, 12.5, 42.5, 47.5, 67.5, 127.5, 132.5, 142.5, 192.5, 217.5))
    expect_equal(map_to_xpos(map, 10),
                 c(5, 15, 45, 55, 75, 135, 145, 155, 205, 230))

})

test_that("map_to_boundaries works", {

    # single chromosome
    expect_equal(map_to_boundaries(list("1"=1:5)), rbind(0,4))

    m <- sort(runif(100, 0, 100))
    expect_equal(map_to_boundaries(list("1"=m)), rbind(0, max(m)-min(m)))

    # two chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85))
    expect_equal(map_to_boundaries(map, 0),
                 cbind(c(0,40), c(40,120)))
    expect_equal(map_to_boundaries(map, 5),
                 cbind(c(0, 45), c(45, 45+80+5)))

    # three chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85),
                "3"=c(90, 100, 150, 175))
    expect_equal(map_to_boundaries(map, 0),
                 cbind(c(0,40), c(40,120), c(120,120+(175-90))))
    expect_equal(map_to_boundaries(map, 5),
                 cbind(c(0,45), c(45,130), c(130,130+(175-90)+5)))
    expect_equal(map_to_boundaries(map, 10),
                 cbind(c(0,50), c(50, 140), c(140, 140+(175-90)+10)))

})

test_that("map_to_index works", {

    # single chromosome
    expect_equal(map_to_index(list("1"=1:5)), list("1"=1:5))

    m <- sort(runif(100, 0, 100))
    expect_equal(map_to_index(list("1"=m)), list("1"=1:100))

    # two chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85))
    expect_equal(map_to_index(map), list("1"=1:3, "2"=4:6))

    # three chromosomes
    map <- list("1"=c(10, 20, 50),
                "2"=c(5,  25, 85),
                "3"=c(90, 100, 150, 175))
    expect_equal(map_to_index(map), list("1"=1:3, "2"=4:6, "3"=7:10))

})

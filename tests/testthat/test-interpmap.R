context("interpolate gen <-> phys maps")

test_that("interpolate_map works", {

    set.seed(31728290)
    oldmap <- sort(runif(20, 0, 100))
    newmap <- sort(runif(20, 0, 200))

    oldpos <- c(-2, seq(0, 100, by=10), 104)

    newpos <- interpolate_map(oldpos, oldmap, newmap)

    # calculate what result should be in R
    expected <- numeric(length(oldpos))
    wh_left <- (oldpos < min(oldmap))
    if(any(wh_left))
        expected[wh_left] <- min(newmap) - (min(oldmap) - oldpos[wh_left])*diff(range(newmap))/diff(range(oldmap))
    wh_right <- (oldpos >= max(oldmap))
    if(any(wh_right))
        expected[wh_right] <- min(newmap) + (oldpos[wh_right] - min(oldmap))*diff(range(newmap))/diff(range(oldmap))
    other <- (!wh_left & !wh_right)
    if(any(other)) {
        index <- sapply(oldpos[other], function(a,b) sum(b <= a), oldmap)
        expected[other] <- newmap[index] + (oldpos[other] - oldmap[index]) *
            (newmap[index+1] - newmap[index]) / (oldmap[index+1] - oldmap[index])
    }

    expect_equal(newpos, expected)


})

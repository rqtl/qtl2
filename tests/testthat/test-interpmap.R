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

test_that("find_intervals works", {
    set.seed(11703070)

    library(qtl)
    data(hyper)
    map <- pull.map(hyper)[[1]]

    pos <- round(runif(1000, 0, 120), 1)
    result <- find_intervals(pos, map)

    expect_true(all(pos[result==-1] < min(map)))
    expect_true(all(pos[result==length(map)-1] >= max(map)))
    middle <- (result >= 0 & result < length(map)-1)
    expect_true(all(pos[middle] >= map[result[middle]+1] &
                    pos[middle] < map[result[middle]+2]))

    # include some positions that are right on the map
    pos <- sample(c(pos, sample(map, 10)))
    result <- find_intervals(pos, map)

    expect_true(all(pos[result==-1] < min(map)))
    expect_true(all(pos[result==length(map)-1] >= max(map)))
    middle <- (result >= 0 & result < length(map)-1)
    expect_true(all(pos[middle] >= map[result[middle]+1] &
                    pos[middle] < map[result[middle]+2]))

    # also check whether is_pos_on_map() works
    on_map <- is_pos_on_map(pos, map, result)
    expect_equal(sum(on_map), 10)
    expect_equal(pos[on_map], map[result[on_map]+1])

})

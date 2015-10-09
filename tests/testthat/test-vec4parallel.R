context("vec4parallel")

test_that("vec4parallel works", {

    # one core
    expect_equivalent(vec4parallel(100, 1), list(1:100))

    # one element
    expect_equivalent(vec4parallel(1, 1), list(1))
    expect_equivalent(vec4parallel(1, 8), list(1))

    # fewer elements than cores
    expect_equivalent(vec4parallel(5, 8), as.list(1:5))
    expect_equivalent(vec4parallel(30, 32), as.list(1:30))

    # equal elements in each core
    expect_equivalent(vec4parallel(10, 5),
                      as.list(as.data.frame(matrix(1:10, ncol=5, byrow=FALSE))))
    expect_equivalent(vec4parallel(144, 8),
                      as.list(as.data.frame(matrix(1:144, ncol=8, byrow=FALSE))))

    # non-equal
    expect_equivalent(vec4parallel(10, 8),
                      list(1:2, 3:4, 5, 6, 7, 8, 9, 10))
    expect_equivalent(vec4parallel(35, 8),
                      list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27, 28:31, 32:35))
    expect_equivalent(vec4parallel(16, 9),
                      list(1:2, 3:4, 5:6, 7:8, 9:10, 11:12, 13:14, 15, 16))

})

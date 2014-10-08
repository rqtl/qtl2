
# this is a helper function for a helper function for read_cross2
context("convert_sexcodes")

test_that("convert_sexcodes works", {

    x <- c(f='female', m='male')
    expected <- c(f=0, m=1)
    expect_equal(convert_sexcodes(x), expected)

    xs <- list(c(f='F', m='M'),
               c(f='Female', m='Male'),
               c(f='FEMALE', m='MALE'),
               c(f='f', m='m'))
    for(x in xs)
        expect_equal(convert_sexcodes(x), expected)

    xs <- list(c(female='female', male='male'),
               c(female='F', male='M'),
               c(female='Female', male='Male'),
               c(female='FEMALE', male='MALE'),
               c(female='f', male='m'))
    expected <- c(female=0, male=1)
    for(x in xs)
        expect_equal(convert_sexcodes(x), expected)

    xs <- list(c("0"='female', "1"='male'),
               c("0"='F', "1"='M'),
               c("0"='Female', "1"='Male'),
               c("0"='FEMALE', "1"='MALE'),
               c("0"='f', "1"='m'))
    expected <- c("0"=0, "1"=1)
    for(x in xs)
        expect_equal(convert_sexcodes(x), expected)

})

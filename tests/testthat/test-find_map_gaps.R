context("find gaps in a map")

test_that("find_map_gaps works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    expect_equal(find_map_gaps(iron$gmap, 40),
                 data.frame(chr=c("1", "4", "5"),
                            left_marker=c("D1Mit80", "D4Mit2", "D5Mit11"),
                            left_index=c(2,1,1),
                            right_marker=c("D1Mit17", "D4Mit352", "D5Mit30"),
                            right_index=c(3,2,2),
                            gap=c(59, 42.7, 44.8),
                            stringsAsFactors=FALSE))

    expect_equal(find_map_gaps(iron$gmap, 38),
                 data.frame(chr=c("1", "4", "5", "12"),
                            left_marker=c("D1Mit80", "D4Mit2", "D5Mit11", "D12Mit88"),
                            left_index=c(2,1,1,1),
                            right_marker=c("D1Mit17", "D4Mit352", "D5Mit30", "D12Mit134"),
                            right_index=c(3,2,2,2),
                            gap=c(59, 42.7, 44.8, 38.2),
                            stringsAsFactors=FALSE))

    expect_equal(find_map_gaps(iron$gmap, 1000),
                 data.frame(chr=character(0),
                            left_marker=character(0),
                            left_index=numeric(0),
                            right_marker=character(0),
                            right_index=numeric(0),
                            gap=numeric(0),
                            stringsAsFactors=FALSE))

})

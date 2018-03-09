context("reduce markers")

test_that("reduce_markers matches qtl::pickMarkerSubset", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    map <- grav2$gmap

    wts <- lapply(vapply(map, length, 1),
                  function(n) runif(n, 1, 5))

    submap <- reduce_markers(map, 1, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 1, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 1, wts[[i]]))
    }

    # repeat with 5 cM
    submap <- reduce_markers(map, 5, wts)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 5, wts[[i]])
        expect_equal(submapi, submap[[i]])

        expect_equal(names(submapi), qtl::pickMarkerSubset(map[[i]], 5, wts[[i]]))
    }

    # without weights, get same result
    set.seed(8764972)
    submap <- reduce_markers(map, 1)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 1)
        expect_equal(submapi, submap[[i]])
    }

    # repeat at 5 cM
    set.seed(8764972)
    submap <- reduce_markers(map, 5)

    # work with each chr one at a time
    #   matches full results?
    #   matches result from R/qtl
    set.seed(8764972)
    for(i in seq(along=map)) {
        submapi <- reduce_markers(map[[i]], 5)
        expect_equal(submapi, submap[[i]])
    }

})


test_that("reduce_markers works", {

    map <- list(c1=c(marker=5),
                c2=c(marker1=5, marker2=7))
    wts <- list(1, c(1,2))
    expect_equal(reduce_markers(map), map)
    expect_equal(reduce_markers(map, 3, wts),
                 list(c1=c(marker=5), c2=c(marker2=7)))

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))

    RNGkind("Mersenne-Twister")
    set.seed(20151205)
    # use proportion typed as weights, but randomize a bit so there are no ties
    prtyped <- lapply(grav2$geno, function(a) colMeans(!is.na(a) & a>0) + runif(ncol(a), 0, 0.1))

    # expected results at min_dist=0.5
    expected <-  structure(list(`1` = structure(c(0, 6.250674, 9.303868, 12.577629,
                     18.39283, 19.648542, 25.23229, 25.882095, 35.669554, 36.631747,
                     44.511901, 50.602339, 53.631519, 57.405778, 62.323561, 67.24062,
                     72.296762, 76.255226, 81.273196, 81.900811, 88.639045, 90.914789,
                     97.121385, 99.711975, 104.80612, 109.519692), .Names = c("PVV4",
                     "AXR-1", "HH.335C-Col/PhyA", "EC.480C", "EC.66C", "GD.86L", "CH.160L-Col",
                     "CC.98L-Col/101C", "AD.121C", "AD.106L-Col", "GB.112L", "GD.97L",
                     "EG.113L/115C", "CD.89C", "BF.206L-Col", "CH.200C", "EC.88C",
                     "GD.160C", "HH.375L", "CH.215L", "BF.116C", "GH.157L-Col", "CC.318C",
                     "CD.173L/175C-Col", "GH.127L-Col/ADH", "HH.360L-Col")), `2` = structure(c(0,
                     6.241796, 11.678237, 12.929907, 13.543347, 14.160665, 14.777988,
                     16.350906, 17.282463, 19.850147, 23.453866, 25.685576, 27.258244,
                     30.511191, 31.128921, 35.085099, 37.651404, 39.551472, 41.451601,
                     42.701883, 47.01718, 47.638238, 49.221621, 52.14781, 55.429739,
                     62.530946), .Names = c("AD.156C", "BF.325L", "GH.580L", "DF.225L",
                     "BH.145C", "FD.226C", "EC.495C-Col", "FD.81L", "BF.105C", "CH.284C",
                     "FD.222L-Col", "CD.245L", "CH.65C", "CH.1500C", "BF.221L", "FD.85C",
                     "GB.150L-Col", "FD.150C", "GD.460L-Col", "Erecta", "BH.195L-Col",
                     "GD.298C", "GH.247L", "BH.120L-Col", "DF.140C", "EG.357C/359L-Col"
                     )), `3` = structure(c(0, 2.613844, 7.032127, 8.291564, 11.480925,
                     15.051063, 16.346984, 18.094189, 18.776988, 21.422012, 22.0854,
                     25.665317, 28.194604, 30.995068, 34.67091, 38.086009, 39.704548,
                     43.418084, 44.350416, 45.169461, 48.653206, 49.640633, 50.30031,
                     52.688081, 55.438502, 58.346106, 59.596912, 66.588704, 68.488354,
                     70.695509, 71.31314, 72.897197, 75.144036, 76.392762), .Names = c("DF.77C",
                     "GB.120C-Col", "GD.248C-Col/249L", "CH.322C", "FD.111L-Col/136C",
                     "CC.266L", "BF.270L-Col/271C", "CC.110L/127C", "EC.58C", "GH.321L/323C-Col",
                     "GH.226C/227L-Col", "EC.83C/84L", "GD.318C/320L", "AD.92L", "HH.410C",
                     "CD.800C", "BF.134C-Col", "EG.122C", "GD.113C-Col", "BH.323C-Col",
                     "DF.303C", "DF.76L", "BF.128C", "HH.117C", "GD.296C-Col", "DF.328C",
                     "FD.98C", "AD.182C", "GD.106C", "HH.171C-Col/173L", "AD.495L-Col",
                     "AD.112L-Col", "BH.109L-Col", "HH.90L-Col")), `4` = structure(c(0,
                     6.388021, 17.889887, 18.5072, 22.755474, 27.797576, 30.442727,
                     34.758259, 35.37558, 39.942822, 40.635774, 45.254346, 46.537827,
                     48.152558, 53.206339, 56.459811, 60.771852, 67.771693, 73.668364
                     ), .Names = c("ANL2", "GH.250C", "GH.165L", "BF.151L", "EC.306L",
                     "BH.92L-Col", "FD.154L", "CH.238C", "CD.84C-Col/85L", "SC5",
                     "DF.108L-Col", "AD.307C", "AD.115L-Col", "FD.167L-Col", "CH.70L/71C-Col",
                     "GH.433L-Col", "GB.490C", "GB.750C", "BH.342C/347L-Col")), `5` = structure(c(0,
                     2.94061, 5.263751, 6.232987, 11.517411, 12.450901, 13.120403,
                     17.456274, 18.082939, 21.390456, 23.320809, 28.452585, 33.209409,
                     34.154297, 38.533295, 41.782964, 45.506147, 51.937046, 53.180697,
                     56.869536, 57.82089, 59.761174, 61.70335, 63.311003, 65.25314,
                     69.302894, 70.2914, 76.061598, 76.680571, 77.933341, 78.551316,
                     79.168573, 80.418873, 81.991829, 82.92342, 85.155018, 87.390432,
                     93.991816, 95.563679, 101.176465, 102.460924, 103.094847, 103.932196,
                     107.488997), .Names = c("FD.207L", "CH.690C", "AD.292L", "BH.144L",
                     "EC.198L-Col", "BH.180C", "BH.325L", "BH.107L-Col", "BF.269C",
                     "AD.114C-Col", "DF.231C", "DF.184L-Col", "GH.473C", "GH.117C",
                     "GH.121L-Col", "DF.154C", "HH.480C", "BH.96L-Col", "CC.188L",
                     "EC.395C", "DF.300C", "GB.235C-Col", "CD.179L", "CH.88L", "CC.216C",
                     "EC.151L", "DFR", "AD.254C", "AD.75C-Col", "GD.350L-Col", "CD.160L",
                     "GB.223C", "DF.450C", "HH.445L-Col", "GD.118C", "CC.262C", "GB.102L-Col/105C",
                     "HH.143C", "BF.168L-Col", "DF.119L", "CH.331L-Col", "GD.222C-Col",
                     "g2368", "EG.205L"))), .Names = c("1", "2", "3", "4", "5"), is_x_chr = structure(c(FALSE,
                     FALSE, FALSE, FALSE, FALSE), .Names = c("1", "2", "3", "4", "5"
                     )))

    attr(expected, "is_x_chr") <- attr(grav2$gmap, "is_x_chr")

    result <- reduce_markers(grav2$gmap, 0.5, prtyped)

    expect_equal(result, expected)

})

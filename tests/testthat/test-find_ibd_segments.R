context("Find IBD segments")

test_that("find_ibd_segments works (regression test)", {

    if(isnt_karl()) skip("this test only run locally")

    recla <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Recla/recla.zip")

    # grab founder genotypes and physical map
    fg <- recla$founder_geno
    pmap <- recla$pmap

    # find shared segments
    segs <- find_ibd_segments(fg, pmap, min_lod=10, error_prob=0.0001)

    expected <-  structure(list(strain1 = c("A", "A", "A", "A", "B", "A", "A",
                                            "A", "A", "A", "B", "A", "A", "D", "C"),
                                strain2 = c("B", "D", "D", "C", "D", "D", "B",
                                            "D", "C", "C", "D", "D", "E", "E", "D"),
                                chr = c("5", "7", "8", "9", "10", "10", "10", "10", "11",
                                        "14", "15", "X", "X", "X", "X"),
                                left = c("backupUNC050583819", "UNC070316647", "backupUNC080772544",
                                         "JAX00688108", "UNC100210555", "JAX00286286",
                                         "JAX00286286", "backupUNC100088273", "UNC110041137",
                                         "UNC141417762", "JAX00059250", "JAX00708714", "UNC200319662",
                                         "UNC200319662", "backupUNC200112735"),
                                right = c("backupUNC050144019", "UNC070374466", "backupUNC080117004",
                                          "UNC090062756", "backupUNC100059055", "backupUNC100059055",
                                          "backupJAX00195786", "JAX00296684", "backupUNC110945346",
                                          "backupUNC140693590", "UNC150504609", "backupUNC200011470",
                                          "backupUNC200388916", "backupUNC200388916", "JAX00719109"),
                                left_pos = c(76.677123, 72.62253, 17.222546, 27.219588, 28.83902,
                                             28.843661, 28.843661, 75.644976, 36.9087, 37.797744,
                                             15.307834, 7.085553, 12.723725, 12.723725, 108.659367),
                                right_pos = c(90.641118, 98.689878, 33.697683, 42.557635, 48.700825,
                                              48.700825, 65.395581, 96.075994, 58.595448, 58.563235,
                                              28.526546, 51.779369, 42.272181, 42.272181, 132.616191),
                                int_length = c(13.963995, 26.067348, 16.475137, 15.338047, 19.861805,
                                               19.857164, 36.55192, 20.431018, 21.686748, 20.765491,
                                               13.218712, 44.693816, 29.548456, 29.548456, 23.956824),
                                n_mar = c(37L, 68L, 39L, 46L, 61L, 60L, 108L, 65L, 62L, 52L, 34L, 96L, 56L,
                                          56L, 53L),
                                n_mismatch = c(0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L,  0L, 0L, 1L, 0L, 0L, 0L),
                                lod = c(11.0161317923257, 12.1048657859049, 10.0495772651155, 10.3967678617075,
                                        11.672980093479, 11.4688763971715, 21.0498587807, 13.6731245696214,
                                        10.2408187280338, 10.7923789614182, 10.606573856941, 15.2505093227274,
                                        10.4956305813686, 10.4956305813686, 11.76581112401)),
                           .Names = c("strain1", "strain2", "chr", "left", "right", "left_pos",
                                      "right_pos", "int_length", "n_mar", "n_mismatch", "lod"),
                           row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                                         "12", "13", "14", "15"), class = "data.frame")

    expect_equal(segs, expected)
})

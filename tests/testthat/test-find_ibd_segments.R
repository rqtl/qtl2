context("Find IBD segments")

test_that("find_ibd_segments works (regression test)", {

    skip_if(isnt_karl(), "this test only run locally")

    recla <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/main/DO_Recla/recla.zip")

    # grab founder genotypes and physical map
    fg <- recla$founder_geno
    pmap <- recla$pmap

    # find shared segments
    segs <- find_ibd_segments(fg, pmap, min_lod=10, error_prob=0.0001)

    expected <-  structure(list(strain1 = c("B", "A", "A", "A", "A", "B", "A", "A",
                                            "A", "A", "A", "B", "A", "A", "D", "C"),
                                strain2 = c("C", "B", "D", "D", "C", "D", "D", "B",
                                            "D", "C", "C", "D", "D", "E", "E", "D"),
                                chr = c("2", "5", "7", "8", "9", "10", "10", "10",
                                        "10", "11", "14", "15", "X", "X", "X", "X"),
                                left_marker = c("JAX00096419", "backupUNC050583819", "UNC070316647",
                                                "JAX00660129", "JAX00688108", "UNC100210555",
                                                "JAX00286286", "JAX00286286", "backupUNC100088273",
                                                "UNC110041137", "UNC141417762", "JAX00059250", "JAX00708714",
                                                "UNC200319662", "UNC200319662", "backupJAX00238910"),
                                right_marker = c("JAX00097346", "backupUNC050144019", "UNC070374466",
                                                 "backupUNC080117004", "UNC090062756", "backupUNC100059055",
                                                 "backupUNC100059055", "backupJAX00195786", "JAX00296684",
                                                 "backupUNC110945346", "backupUNC140693590", "UNC150504609",
                                                 "backupUNC200011470", "backupUNC200388916", "backupUNC200388916",
                                                 "backupUNC200095716"),
                                left_pos = c(80.731655, 76.677123, 72.62253, 8.270413, 27.219588, 28.83902,
                                             28.843661, 28.843661, 75.644976, 36.9087, 37.797744, 15.307834,
                                             7.085553, 12.723725, 12.723725, 41.136597),
                                right_pos = c(93.362628, 90.641118, 98.689878, 33.697683, 42.557635, 48.700825,
                                              48.700825, 65.395581, 96.075994, 58.595448, 58.563235, 28.526546,
                                              51.779369, 42.272181, 42.272181, 66.803155),
                                left_index = c(190L, 178L, 155L, 16L, 65L, 67L, 68L, 68L, 200L, 84L, 75L,
                                               28L, 7L, 24L, 24L, 75L),
                                right_index = c(224L, 214L, 222L, 73L, 110L, 127L, 127L, 175L, 264L, 145L,
                                                126L, 61L, 102L, 79L, 79L, 136L),
                                int_length = c(12.630973, 13.963995, 26.067348, 25.42727, 15.338047, 19.861805,
                                               19.857164, 36.55192, 20.431018, 21.686748, 20.765491, 13.218712,
                                               44.693816, 29.548456, 29.548456, 25.666558),
                                n_mar = c(35L, 37L, 68L, 58L, 46L, 61L, 60L, 108L, 65L, 62L, 52L, 34L, 96L,
                                          56L, 56L, 62L),
                                n_mismatch = c(0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L),
                                lod = c(10.054525618706, 11.0953130382375, 12.162857732679, 14.4366440354073,
                                        10.266444950846, 11.1447388896196, 10.9406351933121, 20.8226367144639,
                                        13.5092786246841, 10.0277601891413, 10.8961283986173, 10.6346025804734,
                                        14.8825542529713, 10.4956305813686, 10.4956305813686, 10.9603560719346)),
                           row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                                         "14", "15", "16"), class = "data.frame")

    expect_equal(segs, expected)
})

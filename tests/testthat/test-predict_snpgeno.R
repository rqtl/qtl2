context("predict_snpgeno")

test_that("predict_snpgeno works", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/DOex/DOex.zip")
    DOex <- read_cross2(file)
    probs <- calc_genoprob(DOex[1:20,"2"], error_prob=0.002)
    m <- maxmarg(probs)

    infg <- predict_snpgeno(DOex, m)
    expect_equal(class(infg), "list")
    expect_equal(length(infg), 1)
    expect_equal(names(infg), "2")
    expect_equal(dim(infg[[1]]), c(20,127))
    expect_true(all(is.na(infg[[1]][,1])))
    expect_equivalent(infg[[1]][9,1:10], c(NA, 2, 1, 1, 3, 2, 2, NA, NA, NA))
    expect_equivalent(infg[[1]][10,101:110], c(2,2,2,2,2,2,2,1,3,2))
    expect_equivalent(infg[[1]][11,51:60], c(NA,NA,NA,NA,1,3,1,2,2,3))

})


test_that("predict_snpgeno works for magic lines", {

    skip_if(isnt_karl(), "this test only run locally")

    # load example data and calculate genotype probabilities
    file <- paste0("https://raw.githubusercontent.com/rqtl/",
                   "qtl2data/main/ArabMAGIC/arabmagic_tair9.zip")
    magic <- read_cross2(file)
    ind <- paste0("MAGIC", ".", 1:20)
    chr <- "2"
    probs <- calc_genoprob(magic[ind,chr], error_prob=0.002)
    m <- maxmarg(probs)

    infg <- predict_snpgeno(magic, m)

    expect_equal(class(infg), "list")
    expect_equal(length(infg), 1)
    expect_equal(names(infg), chr)
    expect_equal(dim(infg[[1]]), c(20,211))
    expect_equivalent(infg[[1]][,1], c(3,1,3,1,NA,3,1,3,1,1,3,3,NA,1,3,NA,1,1,NA,1))
    expect_equivalent(infg[[1]][9,1:10], c(1,1,1,3,3,1,1,3,3,1))
    expect_equivalent(infg[[1]][10,101:110],  c(rep(1,9), 3))
    expect_equivalent(infg[[1]][11,51:60], c(1,1,3,1,1,1,1,1,1,1))

    fg <- magic$founder_geno[[chr]]
    fg[fg==0] <- NA
    infg_hard <- t(sapply(1:20, function(wh_ind) sapply(seq_along(m[[chr]][wh_ind,]), function(i) {
        g <- m[[chr]][wh_ind,i]; ifelse(is.na(g), NA, fg[g,i]) })))
    expect_equivalent(infg[[chr]], infg_hard)

})


test_that("predict_snpgeno works for DOF1 and DOHS", {

    # a small subset of a DOF1 population: 3 individuals, 10 markers on each of 2 chr
    x <- structure(list(crosstype = "dof1", geno = list(`2` = structure(c(3L,
3L, 2L, 3L, 3L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 1L, 1L,
1L, 2L, 2L, 3L, 3L, 3L, 3L, 2L, 2L, 1L, 2L, 2L, 1L), dim = c(3L,
10L), dimnames = list(c("117447", "117443", "119519"), c("UNCHS003833",
"UNCHS003834", "UNC2466989", "JAX00090637", "JAX00482049", "UNCHS003835",
"UNCHS003836", "JAX00482064", "CEAJAX00090652", "JAX00090652"
))), X = structure(c(3L, 3L, 2L, 3L, 3L, 2L, 3L, 3L, 2L, 3L,
3L, 2L, 1L, 1L, 2L, 3L, 3L, 2L, 3L, 3L, 3L, 1L, 1L, 2L, 3L, 3L,
3L, 3L, 3L, 3L), dim = c(3L, 10L), dimnames = list(c("117447",
"117443", "119519"), c("JAX00181486", "JAX00714977", "JAX00181489",
"JAX00714985", "JAX00181502", "JAX00181504", "UNC31049617", "UNCHS048912",
"JAX00715002", "JAX00715002r")))), gmap = list(`2` = c(UNCHS003833 = 0.105908732571738,
UNCHS003834 = 0.10825526663261, UNC2466989 = 0.109617384848426,
JAX00090637 = 0.110597154091381, JAX00482049 = 0.143764599886437,
UNCHS003835 = 0.221542081049763, UNCHS003836 = 0.223175029788022,
JAX00482064 = 0.229358993383802, CEAJAX00090652 = 0.250603258188368,
JAX00090652 = 0.250603258188368), X = c(JAX00181486 = 37.7842866029833,
JAX00714977 = 37.7990420259301, JAX00181489 = 37.8119460221676,
JAX00714985 = 37.8912830161017, JAX00181502 = 37.9925319810839,
JAX00181504 = 37.9928259324011, UNC31049617 = 38.0264384750733,
UNCHS048912 = 38.0297613573523, JAX00715002 = 38.0813853540008,
JAX00715002r = 38.0813853540008)), pmap = list(`2` = c(UNCHS003833 = 3.177758,
UNCHS003834 = 3.181293, UNC2466989 = 3.183345, JAX00090637 = 3.184821,
JAX00482049 = 3.234787, UNCHS003835 = 3.351957, UNCHS003836 = 3.354417,
JAX00482064 = 3.363733, CEAJAX00090652 = 3.395737, JAX00090652 = 3.395737
), X = c(JAX00181486 = 84.865978, JAX00714977 = 84.88179, JAX00181489 = 84.895618,
JAX00714985 = 84.980636, JAX00181502 = 85.089135, JAX00181504 = 85.08945,
UNC31049617 = 85.127743, UNCHS048912 = 85.131962, JAX00715002 = 85.197508,
JAX00715002r = 85.197508)), covar = structure(list(sex = c("M",
"M", "M"), cross_direction = c("(B-B)x(D-D)", "(B-B)x(D-D)",
"(B-B)x(D-D)"), ngen = c("41", "41", "41")), row.names = c("117447",
"117443", "119519"), class = "data.frame"), founder_geno = list(
    `2` = structure(c(1L, 3L, 1L, 3L, 1L, 1L, 3L, 3L, 3L, 3L,
    3L, 3L, 3L, 3L, 1L, 1L, 3L, 3L, 1L, 3L, 1L, 3L, 1L, 3L, 3L,
    1L, 3L, 3L, 1L, 3L, 1L, 3L, 1L, 1L, 3L, 1L, 3L, 3L, 3L, 3L,
    3L, 1L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 1L, 1L, 3L,
    3L, 1L, 3L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 1L,
    3L, 3L, 1L, 1L, 3L, 1L, 3L, 3L, 3L, 1L, 1L, 1L, 1L, 3L, 1L,
    3L, 3L, 3L, 1L, 1L), dim = 9:10, dimnames = list(c("A", "B",
    "C", "D", "E", "F", "G", "H", "I"), c("UNCHS003833", "UNCHS003834",
    "UNC2466989", "JAX00090637", "JAX00482049", "UNCHS003835",
    "UNCHS003836", "JAX00482064", "CEAJAX00090652", "JAX00090652"
    ))), X = structure(c(3L, 1L, 3L, 1L, 1L, 1L, 1L, 1L, 1L,
    3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 1L,
    1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 1L, 3L, 3L, 1L, 1L, 1L,
    1L, 1L, 3L, 3L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 1L, 1L, 3L, 3L,
    3L, 1L, 3L, 1L, 1L, 3L, 3L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L,
    3L, 3L, 1L, 3L, 3L, 3L, 3L, 3L, 1L, 3L, 3L, 3L, 3L, 3L, 3L,
    3L, 3L, 1L, 3L, 3L, 3L), dim = 9:10, dimnames = list(c("A",
    "B", "C", "D", "E", "F", "G", "H", "I"), c("JAX00181486",
    "JAX00714977", "JAX00181489", "JAX00714985", "JAX00181502",
    "JAX00181504", "UNC31049617", "UNCHS048912", "JAX00715002",
    "JAX00715002r")))), is_x_chr = c(`2` = FALSE, X = TRUE),
    is_female = c(`117447` = FALSE, `117443` = FALSE, `119519` = TRUE
    ), cross_info = structure(c(41L, 41L, 41L), dim = c(3L, 1L
    ), dimnames = list(c("117447", "117443", "119519"), "ngen")),
    alleles = c("A", "B", "C", "D", "E", "F", "G", "H")), class = "cross2")

    pr <- calc_genoprob(x, err=0.002)
    set.seed(20221029)
    v <- maxmarg(pr, minprob=0.2)

    snpg <- predict_snpgeno(x, v)

    expect_true(all(snpg$"2" == x$geno$"2"))
    expect_equal(sum(snpg$X != x$geno$X), 2)
    expect_true(all((snpg$X == x$geno$X)[1:2,]))
    expect_equal(snpg$X[3,], colMeans(x$founder_geno$X[c(7,9),]))

})

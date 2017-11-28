context("guess_phase")

test_that("guess_phase works with phase-known data", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
    grav2 <- grav2[1:5,1:2]
    gmap <- insert_pseudomarkers(grav2$gmap, step=5, stepwidth="max")
    pr <- calc_genoprob(grav2, gmap, error_prob=0.002)
    m <- maxmarg(pr)
    expect_equal(guess_phase(grav2, m), m)

})

test_that("guess_phase works with f2 data", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[c(1:5,146:148,211:212), c(19,"X")]
    gmap <- insert_pseudomarkers(iron$gmap, step=5, stepwidth="max")
    pr <- calc_genoprob(iron, gmap, error_prob=0.002)
    set.seed(57407964)
    m <- maxmarg(pr, minprob=0)
    g <- guess_phase(iron, m)

    # autosome
    expect_equal(g[[1]][1,,], cbind(mom=m[[1]][1,]-1, dad=2))
    expect_equal(g[[1]][2,,], cbind(mom=m[[1]][2,]-1, dad=2))
    expect_equal(g[[1]][4,,], cbind(mom=2, dad=m[[1]][2,]-1))
    expect_equal(g[[1]][10,,], cbind(mom=1, dad=m[[1]][10,]))

    # X chromosome
    for(i in 1:5)
        expect_equal(g[[2]][i,,], cbind(mom=m[[2]][i,]-4, dad=NA))
    for(i in 6:8)
        expect_equal(g[[2]][i,,], cbind(mom=m[[2]][i,]-2, dad=2))
    for(i in 9:10)
        expect_equal(g[[2]][i,,], cbind(mom=m[[2]][i,], dad=1))

})



test_that("guess_phase works with DO", {

    set.seed(57407964)
    m <- list("19"=rbind(ind1=c(NA,11,11,11,11,NA,13,13,13,13,13,13,NA,10,10,14,14,NA,10,10),
                         ind2=c(NA,NA, 1,29,32,35,NA,NA,NA,NA, 7,10,10,10,25,19,19,21,19, 7),
                         ind3=c( 2, 2, 2, 8, 8, 8, 8, 8, 8, 8,30,NA, 8,10,10,10, 7, 7, 1, 1)),
              "X"=rbind(ind1=c(43,43,43,43,38,40,40,40,42,42,42,NA,38,38,39,NA,NA,38,38,37),
                        ind2=c(NA,NA,22,22,NA,24, 9,NA,NA,NA,19,19,19,NA,19,19,19,19,19,16),
                        ind3=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)))

    for(i in 1:2)
        colnames(m[[i]]) <- paste0("marker", 1:20)

    # fake cross2 object
    x <- list(crosstype="do",
              geno=list("19"=rbind(ind1=0, ind2=0, ind3=0), X=rbind(ind1=0, ind2=0, ind3=0)),
              is_female=c(ind1=FALSE, ind2=TRUE, ind3=TRUE),
              is_x_chr=c("19"=FALSE,X=TRUE))
    class(x) <- c("cross2", "list")

    g <- guess_phase(x, m)

    expected <- list("19"=array(dim=c(3,20,2)),
                     X=array(dim=c(3,20,2)))
    for(i in 1:2)
        dimnames(expected[[i]]) <- list(c("ind1", "ind2","ind3"),
                                        colnames(m[[1]]),
                                        c("mom", "dad"))


    expected[[1]][1,,1] <- c(NA,1,1,1,1,NA,3,3,3,3,3,3,NA,4,4,4,4,NA,4,4)
    expected[[1]][1,,2] <- c(NA,5,5,5,5,NA,5,5,5,5,5,5,NA,4,4,5,5,NA,4,4)
    expected[[1]][2,,1] <- c(NA,NA,1,1,4,7,NA,NA,NA,NA,1,4,4,4,4,4,4,6,4,4)
    expected[[1]][2,,2] <- c(NA,NA,1,8,8,8,NA,NA,NA,NA,4,4,4,4,7,6,6,6,6,1)
    expected[[1]][3,,1] <- c(1,1,1,4,4,4,4,4,4,4,8,NA,4,4,4,4,1,1,1,1)
    expected[[1]][3,,2] <- c(2,2,2,2,2,2,2,2,2,2,2,NA,2,4,4,4,4,4,1,1)


    expected[[2]][1,,1] <- m[[2]][1,]-36
    expected[[2]][2,,1] <- c(NA,NA,7,7,NA,7,4,NA,NA,NA,4,4,4,NA,4,4,4,4,4,1)
    expected[[2]][2,,2] <- c(NA,NA,1,1,NA,3,3,NA,NA,NA,6,6,6,NA,6,6,6,6,6,6)

    expect_equal(g, expected)

    # deterministic version
    set.seed(57407965) # the seed shouldn't matter here
    g <- guess_phase(x, m, deterministic=TRUE)

    # expected when deterministic: just one change
    expected[[2]][2,,1:2] <- expected[[2]][2,,2:1]

    expect_equal(g, expected)

})

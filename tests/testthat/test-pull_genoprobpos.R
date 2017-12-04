context("pull genotype probabilities for a specified position")

test_that("pull_genoprobpos works", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    iron <- iron[,c(2,8,"X")]
    gmap <- insert_pseudomarkers(iron$gmap, step=1)
    pr <- calc_genoprob(iron, gmap, error_prob=0.002)

    pmar <- find_marker(gmap, 8, 40)
    pr_8_40 <- pull_genoprobpos(pr, pmar)
    expect_equal(pr_8_40, pr[["8"]][,,"c8.loc40"])

    expect_equal(pull_genoprobpos(pr["275",], pmar),
                 pr_8_40["275",,drop=FALSE])

    expect_warning(pull_genoprobpos(pr, c("c8.loc40", find_marker(gmap, 8, 70))))

    expect_error(pull_genoprobpos(pr, "karl"))

})

context("is_phase_known")

test_that("is_phase_known works for grav2 and iron", {

    grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    expect_true(is_phase_known(grav2))

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    expect_false(is_phase_known(iron))

})


test_that(".is_phase_known works for all cross types", {

    for(crosstype in c("bc", "haploid", "dh", "riself", "risib",
                       "riself4", "riself8", "riself16", "risib4", "risib8",
                       "magic19", "dh6", "dof1")) {
        expect_true(.is_phase_known(crosstype))
    }

    for(crosstype in c("f2", "ail", "ail3", "do", "hs")) {
        expect_false(.is_phase_known(crosstype))
    }

})

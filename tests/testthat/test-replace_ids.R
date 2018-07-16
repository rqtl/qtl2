context("replace individual IDs")

test_that("check_new_ids() works", {

    old_ids <- as.character(1:200)
    new_ids <- setNames(sprintf("mouse%03d", as.numeric(old_ids)), old_ids)

    # same ids in same order
    expect_equal(check_new_ids(new_ids, old_ids), new_ids)

    # same ids, but shuffled
    shuffled <- sample(new_ids)
    expect_equal(check_new_ids(shuffled, old_ids), shuffled)

    # error if duplicate IDs
    dup_new <- new_ids
    dup_new[5] <- dup_new[20]
    expect_error(check_new_ids(dup_new, old_ids))

    # error if duplicate names in the IDs
    dup_old <- new_ids
    names(dup_old)[5] <- names(dup_old)[20]
    expect_error(check_new_ids(dup_old, old_ids))

    # warning if extra IDs
    new_extra <- c(new_ids, "201"="mouse201", "202"="mouse202")
    expect_warning( expect_equal( check_new_ids(new_extra, old_ids), new_ids) )
    o <- sample(length(new_extra))
    expect_warning( expect_equal( check_new_ids(new_extra[o], old_ids),
                                  new_extra[o[o<=length(new_ids)]] ) )

    # warning if not all IDs are there
    new_missing <- sample(new_ids, length(new_ids)-5)
    expect_warning( expect_equal( check_new_ids(new_missing, old_ids),
                                  new_missing[names(new_missing) %in% old_ids] ))


})


test_that("replace_ids() works for a cross2 object", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    ids <- ind_ids(iron)
    new_ids <- setNames(paste0("mouse", ids), ids)
    change_back <- setNames(ids, paste0("mouse", ids))

    # same ids, old and new (replace then replace back)
    expect_equal( replace_ids(replace_ids(iron, new_ids), change_back), iron)

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})


test_that("replace_ids() works for calc_genoprob output", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    ids <- ind_ids(iron)
    new_ids <- setNames(paste0("mouse", ids), ids)
    change_back <- setNames(ids, paste0("mouse", ids))

    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    pr <- calc_genoprob(iron, map)

    # same ids, old and new (changed back)
    expect_equal( replace_ids(replace_ids(pr, new_ids), change_back), pr)

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})


test_that("replace_ids() works for viterbi output", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    ids <- ind_ids(iron)
    new_ids <- setNames(paste0("mouse", ids), ids)
    change_back <- setNames(ids, paste0("mouse", ids))

    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    v <- viterbi(iron, map)

    # same ids, old and new (changed back)
    expect_equal( replace_ids(replace_ids(v, new_ids), change_back), v)

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})


test_that("replace_ids() works for sim_geno output", {

    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
    ids <- ind_ids(iron)
    new_ids <- setNames(paste0("mouse", ids), ids)
    change_back <- setNames(ids, paste0("mouse", ids))

    map <- insert_pseudomarkers(iron$gmap, step=2.5)
    d <- sim_geno(iron, map, n_draws=8)

    # same ids, old and new (changed back)
    expect_equal( replace_ids(replace_ids(d, new_ids), change_back), d)

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})

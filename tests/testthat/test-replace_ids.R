context("replace individual IDs")

test_that("check_new_ids() works", {

    old_ids <- as.character(1:200)
    new_ids <- setNames(sprintf("mouse%03d", as.numeric(old_ids)), old_ids)

    # same ids, old and new
    expect_equal(check_new_ids(new_ids, old_ids), new_ids)

    # error if duplicate IDs
    dup_new <- new_ids
    dup_new[5] <- dup_new[20]
    expect_error(check_new_ids(dup_new, old_ids))

    # error if duplicate names in the IDs
    dup_old <- new_ids
    names(dup_old)[5] <- names(dup_old)[20]
    expect_error(check_new_ids(dup_old, old_ids))

    # warning if extra IDs

    # warning if not all IDs are there

    # simple replacement, everything in order

    # simple replacement, but shuffled



})


test_that("replace_ids() works for a cross2 object", {

    # same ids, old and new

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})


test_that("replace_ids() works for calc_genoprob output", {

    # same ids, old and new

    # simple replacement, everything in order

    # simple replacement, but shuffled

    # simple replacement, with some extras plus shuffled

    # missing some individuals


})

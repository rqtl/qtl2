#' Replace individual IDs
#'
#' Replace the individual IDs in an object with new ones
#'
#' @param x Object whose IDs will be replaced
#'
#' @param ids Vector of character strings with the new individual IDs, with the names being the original IDs.
#'
#' @return The input `x` object, but with individual IDs replaced.
#'
#' @importFrom stats setNames
#' @export
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' ids <- as.numeric(ind_ids(iron))
#'
#' # replace the numeric IDs with IDs like "mouse003"
#' new_ids <- setNames( sprintf("mouse%03d", as.numeric(ids)), ids)
#'
#' iron <- replace_ids(iron, new_ids)
replace_ids <- function(x, ids) UseMethod("replace_ids")


# check that the new ids vector is okay
check_new_ids <-
    function(new_ids, old_ids)
{
    if(length(new_ids) != length(unique(new_ids)))
        stop("The new IDs are not unique")

    new_ids_names <- names(new_ids)
    if(is.null(new_ids_names)) stop("The vector of new ids needs names, indicating the corresponding old IDs")
    if(length(new_ids_names) != length(unique(new_ids_names)))
        stop("The names of the new IDs are not unique")

    if(length(old_ids) != length(unique(old_ids)))
        stop("The old IDs are not unique")

    if(!all(new_ids_names %in% old_ids)) {
        warning("The new IDs include some extras")
        new_ids <- new_ids[new_ids_names %in% old_ids]
        new_ids_names <- names(new_ids)
    }

    if(!all(old_ids %in% new_ids_names)) {
        warning("Some old IDs are missing from the new IDs object")
    }

    new_ids
}


#' @describeIn replace_ids Replace IDs in a `"cross2"` object
#' @export
replace_ids.cross2 <-
    function(x, ids)
{
    ids <- check_new_ids(ids, ind_ids(x))
    ids_names <- names(ids)
    names(ids) <- NULL

    # rownames for each chr in geno
    for(i in seq_along(x$geno)) {
        m <- match(rownames(x$geno[[i]]), ids_names)
        if(any(is.na(m))) x$geno[[i]] <- x$geno[[i]][!is.na(m),,drop=FALSE]
        rownames(x$geno[[i]]) <- ids[m[!is.na(m)]]
    }

    # rownames in pheno, covar, cross_info
    for(obj in c("pheno", "covar", "cross_info")) {
        if(!(obj %in% names(x))) next

        m <- match(rownames(x[[obj]]), ids_names)
        if(any(is.na(m))) x[[obj]] <- x[[obj]][!is.na(m),,drop=FALSE]
        rownames(x[[obj]]) <- ids[m[!is.na(m)]]
    }

    # names in is_female
    for(obj in c("is_female")) {
        if(!obj %in% names(x)) next

        m <- match(names(x[[obj]]), ids_names)
        if(any(is.na(m))) x[[obj]] <- x[[obj]][!is.na(m)]
        names(x[[obj]]) <- ids[m[!is.na(m)]]
    }

    x
}




#' @describeIn replace_ids Replace IDs in output from [calc_genoprob()]
#' @export
replace_ids.calc_genoprob <-
    function(x, ids)
{
    ids <- check_new_ids(ids, rownames(x[[1]]))
    ids_names <- names(ids)
    names(ids) <- NULL

    # replace the row names for each chromosome; possibly subsetting things
    for(i in seq_along(x)) {
        m <- match(rownames(x[[i]]), ids_names)
        x[[i]] <- x[[i]][!is.na(m),,,drop=FALSE]
        rownames(x[[i]]) <- ids[m[!is.na(m)]]
    }

    x
}

#' @describeIn replace_ids Replace IDs in output from [viterbi()]
#' @export
replace_ids.viterbi <-
    function(x, ids)
{
    ids <- check_new_ids(ids, rownames(x[[1]]))
    ids_names <- names(ids)
    names(ids) <- NULL

    # replace the row names for each chromosome; possibly subsetting things
    for(i in seq_along(x)) {
        m <- match(rownames(x[[i]]), ids_names)
        x[[i]] <- x[[i]][!is.na(m),,drop=FALSE]
        rownames(x[[i]]) <- ids[m[!is.na(m)]]
    }

    x

}

#' @describeIn replace_ids Replace IDs in output from [sim_geno()]
#' @export
replace_ids.sim_geno <- function(x, ids) replace_ids.calc_genoprob(x, ids)

#' @describeIn replace_ids Replace IDs in a matrix
#' @export
replace_ids.matrix <-
    function(x, ids)
{
    ids <- check_new_ids(ids, rownames(x))
    ids_names <- names(ids)
    names(ids) <- NULL

    m <- match(rownames(x), ids_names)
    x <- x[!is.na(m),,drop=FALSE]
    rownames(x) <- ids[m[!is.na(m)]]

    x
}

#' @describeIn replace_ids Replace IDs in a data frame
#' @export
replace_ids.data.frame <- replace_ids.matrix

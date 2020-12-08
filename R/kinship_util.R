# check that kinship matrices are square
# and have same row and column IDs
#
# returns vector of IDs
check_kinship <-
    function(kinship, n_chr)
{
    if(is_kinship_decomposed(kinship)) {
        if(is_kinship_list(kinship)) {
            if(length(kinship) != n_chr)
                stop("length(kinship) != no. chromosomes (", n_chr, ")")
            return(rownames(kinship[[1]]$vectors))
        }
        else {
            return(rownames(kinship$vectors))
        }
    }

    if(!is_kinship_list(kinship)) { # one kinship matrix
        stopifnot(nrow(kinship) == ncol(kinship))
        stopifnot( all(rownames(kinship) == colnames(kinship)) )
        return(rownames(kinship))
    } else {
        if(length(kinship) != n_chr)
            stop("length(kinship) != no. chromosomes (", n_chr, ")")

        kinship_square <- vapply(kinship, function(mat) nrow(mat) == ncol(mat), TRUE)
        stopifnot( all(kinship_square) )

        kinship_sameIDs <- vapply(kinship, function(mat) (nrow(mat) == nrow(kinship[[1]])) &&
                                  all((rownames(mat) == rownames(kinship[[1]])) &
                                      (colnames(mat) == colnames(kinship[[1]])) &
                                      (rownames(mat) == colnames(kinship[[1]]))), TRUE)
        if(!all(kinship_sameIDs))
            stop("All kinship matrices should be the same size ",
                 "and have the same row and column names")

        return(rownames(kinship[[1]]))
    }
}

# multiply kinship matrix by 2
# see Almasy & Blangero (1998) https://doi.org/10.1086/301844
#
# This can also handle the case of "loco", and of having eigen decomposition pre-computed
double_kinship <-
    function(kinship)
{
    if(is.null(kinship)) return(NULL)

    if(is_kinship_decomposed(kinship)) { # already did decomposition
        if(is_kinship_list(kinship)) { # list of decomposed kinship matrices
            kinship <- lapply(kinship, function(a) { a$values <- 2*a$values; a })
        }
        else { # one decomposed kinship matrix
            kinship$values <- 2*kinship$values
        }
    }
    else {
        if(is.list(kinship))
            kinship <- lapply(kinship, function(a) a*2)
        else kinship <- 2*kinship
    }

    kinship
}


# check if alread decomposed
is_kinship_decomposed <-
    function(kinship)
{
    decomp <- attr(kinship, "eigen_decomp")

    (!is.null(decomp) && decomp) || # should have attribute
        (length(kinship)==2 && all(names(kinship) == c("values", "vectors"))) || # single-chr case missing attribute
        (is.list(kinship) && all(vapply(kinship, length, 1)==2) && all(vapply(kinship, function(a) all(names(a)==c("values", "vectors")), TRUE))) # multi-chr case
}

# is kinship a list with (potentially) multiple chromosomes
is_kinship_list <-
    function(kinship)
{
    if(is_kinship_decomposed(kinship)) {
        if(length(kinship)==2 && all(names(kinship) == c("values", "vectors"))) { # just one chromosome
            return(FALSE)
        }
        else return(TRUE)
    }
    else {
        return(is.list(kinship))
    }
}

# check that kinship concerns one chromosome
# and remove outer list if necessary
check_kinship_onechr <-
    function(kinship)
{
    if(is_kinship_list(kinship)) {
        if(length(kinship) > 1)
            stop("kinship should concern just one chromosome")
        decomp <- attr(kinship, "eigen_decomp") ## preserve attribute
        kinship <- kinship[[1]]
        attr(kinship, "eigen_decomp") <- decomp
    }

    kinship
}

# multiply kinship by weights (from and back)
# assuming weights are really square-root weights
#' @importFrom stats setNames
weight_kinship <-
    function(kinship, weights=NULL, tol=1e-8)
{
    # if null weights are all 1's, just skip it
    if(is.null(weights) || max(abs(weights-1)) < tol) return(kinship)

    if(is_kinship_list(kinship)) {
        for(i in seq_along(kinship)) {
            kinship[[i]] <- weight_kinship(kinship[[i]], weights)
        }
        return(kinship)
    }

    # if kinship was decomposed, expand it and then decompose it later
    do_decomp <- FALSE
    if(is_kinship_decomposed(kinship)) {
        do_decomp <- TRUE
        # expand out the decomposition
        kinship <- t(kinship$vectors) %*% diag(kinship$values) %*% kinship$vectors
    }

    # line them up
    ind2keep <- get_common_ids(setNames(rownames(kinship), NULL), setNames(names(weights), NULL))
    weights <- weights[ind2keep]
    kinship <- kinship[ind2keep, ind2keep, drop=FALSE]

    # multiply kinship matrix by weights
    kinship <- kinship * outer(weights, weights, "*")

    # decompose kinship again
    if(do_decomp) kinship <- decomp_kinship(kinship)

    kinship
}

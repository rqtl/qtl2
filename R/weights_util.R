# utilities for weights

# are the weights NULL or all 0's?
is_null_weights <-
    function(weights, tol=1e-12)
{
    if(is.null(weights) || max(abs(weights-1)) < tol)
        return(TRUE)

    FALSE
}


# check vector of weights & take square-roots
#
# - all positive?
#
# - if no missing values and all close to 1,
#   just use NULL rather than the weights

sqrt_weights <-
    function(weights, tol=1e-12)
{
    if(is.null(weights)) return(weights)

    if(any(weights <= 0))
        stop("weights must all be positive")

    if(all(!is.na(weights) & abs(weights - 1)<tol))
        return(NULL)

    weights <- stats::setNames( as.numeric(weights), names(weights) )

    return(sqrt(weights))
}


# multiply a vector by a set of weights
weight_vector <-
    function(vec, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(vec)) return(vec)

    # align and multiply
    id <- get_common_ids(names(vec), names(weights))
    vec[id] * weights[id]
}


# multiply a matrix by a set of weights
weight_matrix <-
    function(mat, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(mat)) return(mat)

    # force the input mat to be a matrix
    if(!is.matrix(mat)) mat <- as.matrix(mat)

    # align and multiply
    id <- get_common_ids(rownames(mat), names(weights))
    mat[id,,drop=FALSE] * weights[id]
}

# multiply an array by a set of weights
weight_array <-
    function(arr, weights, tol=1e-12)
{
    if(is_null_weights(weights, tol) || is.null(arr)) return(arr)

    # force the input mat to be a matrix
    if(!is.array(arr) || length(dim(arr)) != 3)
        stop("arr should be a 3-dimensional array")

    # align and multiply
    id <- get_common_ids(rownames(arr), names(weights))
    arr[id,,,drop=FALSE] * weights[id]
}

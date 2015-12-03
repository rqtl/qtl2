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

    return(sqrt(weights))
}

#' Scale kinship matrix
#'
#' Scale kinship matrix to be like a correlation matrix.
#'
#' @param kinship A kinship matrix, or a list of such in the case of
#'     the "leave one chromosome out" method, as calculated by
#'     [calc_kinship()].
#'
#' @return A matrix or list of matrices, as with the input, but with
#'     the matrices scaled to be like correlation matrices.
#'
#' @details We take
#' \ifelse{html}{\out{<em>c<sub>ij</sub> = k<sub>ij</sub> / &radic;(k<sub>ii</sub> k<sub>jj</sub>)</em>}}{\eqn{c_{ij} = k_{ij} / \sqrt{k_{ii} k_{jj}}}{c_ij = k_ij / sqrt(k_ii k_jj)}}
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' map <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, map, error_prob=0.002)
#' K <- calc_kinship(probs)
#' Ka <- scale_kinship(K)

scale_kinship <-
    function(kinship)
{
    if(is_kinship_decomposed(kinship))
        stop("scale_kinship cannot be applied, because kinship has been decomposed")

    k_attr <- attributes(kinship)

    if(is.list(kinship))
        result <- lapply(kinship, scale_kinship)
    else {
        d <- diag(kinship)
        if(any(d <= 0)) stop("Some diagonal elements are <= 0")
        d <- sqrt(d)
        result <- t(t(kinship/d) / d)
    }

    for(a in names(k_attr)) attr(result, a) <- k_attr[[a]]

    result
}

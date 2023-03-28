#' Unsmooth genetic map
#'
#' Performs the reverse operation of [smooth_gmap()], in case one wants to go back
#' to the original genetic map.
#'
#' @param gmap Genetic map, as a list of numeric vectors; each vector gives marker
#' positions for a single chromosome.
#' @param pmap Physical map, as a list of numeric vectors; each vector gives marker
#' positions for a single chromosome, with the same chromosomes and markers as `gmap`.
#' @param alpha Proportion of mixture to take from constant recombination.
#'
#' @return A genetic map like the input `gmap`, but with the reverse
#' operation of [smooth_gmap()] applied, provided that exactly the
#' same physical map and `alpha` are used.
#'
#' @details An interval of genetic length \eqn{d_g}{dg} and physical
#'     length \eqn{d_p}{dp} is changed to have length
#'     \eqn{(d_g - \alpha d_p r)/(1-\alpha)}{(dg + alpha*dp*r)/(1-alpha)}
#'     where \eqn{r = L_g / L_p}{r = Lg/Lp} is the chromosome-specific
#'     recombination rate.
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' gmap_adj <- smooth_gmap(iron$gmap, iron$pmap)
#' gmap_back <- unsmooth_gmap(gmap_adj, iron$pmap)
#'
#' @export
#' @keywords utilities
#' @seealso [smooth_gmap()]

unsmooth_gmap <-
    function(gmap, pmap, alpha=0.02)
    {
        stopifnot(is.list(gmap), is.list(pmap))
        stopifnot(length(gmap)==length(pmap), all(names(gmap) == names(pmap)))
        stopifnot(length(alpha)==1, alpha >= 0, alpha <= 1)
        for(chr in seq(along=gmap)) {
            g <- gmap[[chr]]
            p <- pmap[[chr]]
            stopifnot(length(g) == length(p), all(names(g) == names(p)) )
            gd <- diff(g)
            pd <- diff(p)
            rr <- diff(range(g)) / diff(range(p))
            gd <- (gd- alpha*pd*rr)/(1-alpha)
            gmap[[chr]] <- setNames(cumsum(c(g[1], gd)), names(g))
        }
        gmap
    }

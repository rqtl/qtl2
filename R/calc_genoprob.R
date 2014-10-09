# calc_genoprob
#' Calculate conditional genotype probabilities
#'
#' Uses a hidden Markov model to calculate the probabilities of the
#' true underlying genotypes given the observed multipoint marker
#' data, with possible allowance for genotyping errors.
#'
#' @param cross
#' @param step Distance between pseudomarkers and markers; if
#' \code{step=0} no pseudomarkers are inserted.
#' @param off_end Distance beyond terminal markers in which to insert
#' pseudomarkers.
#' @param stepwidth Indicates whether to use a fixed grid
#' (\code{stepwidth="fixed"}) or to use the maximal distance between
#' pseudomarkers to ensure that no two adjacent markers/pseudomarkers
#' are more than \code{step} apart.
#' @param pseudomarker_map A map of pseudomarker locations; if provided the
#' \code{step}, \code{off_end}, and \code{stepwidth} arguments are
#' ignored.
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#'
#' @return A list of three-dimensional arrays of probabilities,
#' individuals x positions x genotypes
#'
#' @details
#'   Let \eqn{O_k}{O[k]} denote the observed marker genotype at position
#'  \eqn{k}, and \eqn{g_k}{g[k]} denote the corresponding true underlying
#'  genotype.
#'
#'  We use the forward-backward equations to calculate
#'  \eqn{\alpha_{kv} = \log Pr(O_1, \ldots, O_k, g_k = v)}{%
#'    a[k][v] = log Pr(O[1], \ldots, O[k], g[k] = v)}
#'  and
#'  \eqn{\beta_{kv} = \log Pr(O_{k+1}, \ldots, O_n | g_k = v)}{%
#'    b[k][v] = log Pr(O[k+1], \ldots, O[n] | g[k] = v)}
#'
#'  We then obtain
#'  \eqn{Pr(g_k | O_1, \ldots, O_n) = \exp(\alpha_{kv} + \beta_{kv}) / s}{%
#'    Pr(g[k] | O[1], \ldots, O[n] = exp(a[k][v] + b[k][v]) / s}
#'  where
#'  \eqn{s = \sum_v \exp(\alpha_{kv} + \beta_{kv})}{%
#'    s = sum_v exp(a[k][v] + b[k][v])}
#'
#' @export
#' @keywords utilities
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

calc_genoprob <-
function(cross, step=0, off_end=0, stepwidth=c("fixed", "max"), pseudomarker_map,
         error_prob=1e-4, map_function=c("haldane", "kosambi", "c-f", "morgan"))
{
    # check inputs
    if(class(cross) != "cross2")
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)
    stepwidth <- match.arg(stepwidth)

    # construct map at which to do the calculations
    if(missing(pseudomarker_map))
        pseudomarker_map <- NULL
    # tolerance for matching marker and pseudomarker positions
    tol <- ifelse(step==0 || step>1, 0.01, step/100)
    # create the combined marker/pseudomarker map
    map <- insert_pseudomarkers(cross$gmap, step, off_end, stepwidth,
                                pseudomarker_map, tol)

    probs <- vector("list", length(map))
    rf <- lapply(map, function(m) mf(diff(m), map_function))
    cross_info <- t(cross$cross_info)

    for(i in seq(along=map)) {
        probs[[i]] <- .calc_genoprob(cross$crosstype, t(cross$geno[[i]]),
                                     cross$is_x_chr[i], cross$is_female,
                                     cross_info, rf[[i]], attr(map[[i]], "index"),
                                     error_prob) %>% aperm(c(2,3,1))

        dimnames(probs[[i]]) <- list(rownames(cross$geno),
                                     names(map[[i]]),
                                     NULL) # FIX_ME: need genotype names in here
    }

    names(probs) <- names(cross$gmap)
    attr(probs, "map") <- map

    probs
}

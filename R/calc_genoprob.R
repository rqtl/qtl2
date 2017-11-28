# calc_genoprob
#' Calculate conditional genotype probabilities
#'
#' Uses a hidden Markov model to calculate the probabilities of the
#' true underlying genotypes given the observed multipoint marker
#' data, with possible allowance for genotyping errors.
#'
#' @md
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param map Genetic map of markers. May include pseudomarker
#' locations (that is, locations that are not within the marker
#' genotype data). If NULL, the genetic map in `cross` is used.
#' @param error_prob Assumed genotyping error probability
#' @param map_function Character string indicating the map function to
#' use to convert genetic distances to recombination fractions.
#' @param lowmem If `FALSE`, split individuals into groups with
#' common sex and crossinfo and then precalculate the transition
#' matrices for a chromosome; potentially a lot faster but using more
#' memory.
#' @param quiet If `FALSE`, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A list of three-dimensional arrays of probabilities,
#'     individuals x genotypes x positions. (Note that the arrangement is
#'     different from R/qtl.) Also contains four attributes:
#' * `crosstype` - The cross type of the input `cross`.
#' * `is_x_chr` - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input `cross`.
#' * `alleles` - Vector of allele codes, from input
#'     `cross`.
#' * `alleleprobs` - Logical value (`FALSE`) that
#'     indicates whether the probabilities are compressed to allele
#'     probabilities, as from [genoprob_to_alleleprob()].
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
#' @seealso [insert_pseudomarkers()]
#'
#' @examples
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2"))
#' gmap_w_pmar <- insert_pseudomarkers(grav2$gmap, step=1)
#' probs <- calc_genoprob(grav2, gmap_w_pmar, error_prob=0.002)

calc_genoprob <-
function(cross, map=NULL, error_prob=1e-4,
         map_function=c("haldane", "kosambi", "c-f", "morgan"),
         lowmem=FALSE, quiet=TRUE, cores=1)
{
    # check inputs
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)

    if(!lowmem) # use other version
        return(calc_genoprob2(cross=cross, map=map,
                              error_prob=error_prob, map_function=map_function,
                              quiet=quiet, cores=cores))

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # pseudomarker map
    if(is.null(map))
        map <- insert_pseudomarkers(cross$gmap)
    # possibly subset the map
    if(length(map) != length(cross$geno) || !all(names(map) == names(cross$geno))) {
        chr <- names(cross$geno)
        if(!all(chr %in% names(map)))
            stop("map doesn't contain all of the necessary chromosomes")
        map <- map[chr]
    }
    # calculate marker index object
    index <- create_marker_index(lapply(cross$geno, colnames), map)

    probs <- vector("list", length(map))
    rf <- map2rf(map, map_function)

    # deal with missing information
    ind <- rownames(cross$geno[[1]])
    chrnames <- names(cross$geno)
    is_x_chr <- handle_null_isxchr(cross$is_x_chr, chrnames)
    cross$is_female <- handle_null_isfemale(cross$is_female, ind)
    cross$cross_info <- handle_null_isfemale(cross$cross_info, ind)

    founder_geno <- cross$founder_geno
    if(is.null(founder_geno))
        founder_geno <- create_empty_founder_geno(cross$geno)

    by_group_func <- function(i) {
        pr <- .calc_genoprob(cross$crosstype, t(cross$geno[[chr]][group[[i]],,drop=FALSE]),
                             founder_geno[[chr]], cross$is_x_chr[chr], cross$is_female[group[[i]]],
                             t(cross$cross_info[group[[i]],,drop=FALSE]), rf[[chr]], index[[chr]],
                             error_prob)
        aperm(pr, c(2,1,3))
    }

    group <- parallel::splitIndices(nrow(cross$geno[[1]]), n_cores(cores))
    groupindex <- seq(along=group)

    probs <- vector("list", length(cross$geno))
    for(chr in seq(along=cross$geno)) {
        if(!quiet) message("Chr ", names(cross$geno)[chr])

        if(n_cores(cores) == 1) { # no parallel processing
            # calculations in one group
            probs[[chr]] <- by_group_func(1)
        }
        else {
            # calculations in parallel
            temp <- cluster_lapply(cores, groupindex, by_group_func)

            # paste them back together
            d <- vapply(temp, dim, rep(0,3))
            nr <- sum(d[1,])
            probs[[chr]] <- array(dim=c(nr, d[2,1], d[3,1]))
            for(i in groupindex)
                probs[[chr]][group[[i]],,] <- temp[[i]]
        }

        # genotype names
        alleles <- cross$alleles
        if(is.null(alleles)) # no alleles saved; use caps
            alleles <- LETTERS[1:ncol(probs[[chr]])]
        gnames <- geno_names(cross$crosstype,
                             alleles,
                             cross$is_x_chr[chr])
        if(length(gnames) != ncol(probs[[chr]])) {
            warning("genotype names has length (", length(gnames),
                    ") != genoprob dim (", ncol(probs[[chr]]), ")")
            gnames <- NULL
        }

        dimnames(probs[[chr]]) <- list(rownames(cross$geno[[chr]]),
                                       gnames,
                                       names(map[[chr]]))

    }

    names(probs) <- names(cross$gmap)

    attr(probs, "crosstype") <- cross$crosstype
    attr(probs, "is_x_chr") <- cross$is_x_chr
    attr(probs, "alleles") <- cross$alleles
    attr(probs, "alleleprobs") <- FALSE

    class(probs) <- c("calc_genoprob", "list")

    probs
}

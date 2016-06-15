# calc_genoprob
#' Calculate conditional genotype probabilities
#'
#' Uses a hidden Markov model to calculate the probabilities of the
#' true underlying genotypes given the observed multipoint marker
#' data, with possible allowance for genotyping errors.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
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
#' @param lowmem If \code{FALSE}, split individuals into groups with
#' common sex and crossinfo and then precalculate the transition
#' matrices for a chromosome; potentially a lot faster but using more
#' memory.
#' @param quiet If \code{FALSE}, print progress messages.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return A list containing the following
#' \itemize{
#' \item \code{probs} - A list of three-dimensional arrays of probabilities,
#'     individuals x genotypes x positions. (Note that the arrangement is
#'     different from R/qtl.)
#' \item \code{map} - The genetic map as a list of vectors of marker positions.
#' \item \code{grid} - A list of logical vectors, indicating which
#'     positions correspond to a grid of markers/pseudomarkers. (may be
#'     absent)
#' \item \code{crosstype} - The cross type of the input \code{cross}.
#' \item \code{is_x_chr} - Logical vector indicating whether chromosomes
#'     are to be treated as the X chromosome or not, from input \code{cross}.
#' \item \code{is_female} - Vector of indicators of which individuals are female, from input
#'     \code{cross}.
#' \item \code{cross_info} - Matrix of cross information for the
#'     individuals, from input \code{cross}.
#' \item \code{alleles} - Vector of allele codes, from input
#'     \code{cross}.
#' \item \code{alleleprobs} - Logical value (\code{FALSE}) that
#'     indicates whether the probabilities are compressed to allele
#'     probabilities, as from \code{\link{genoprob_to_alleleprob}}.
#' \item \code{step} - the value of the \code{step} argument.
#' \item \code{off_end} - the value of the \code{off_end} argument.
#' \item \code{stepwidth} - the value of the \code{stepwidth} argument.
#' \item \code{error_prob} - the value of the \code{error_prob} argument.
#' \item \code{map_function} - the value of the \code{map_function} argument.
#' }
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
#' grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
#' probs <- calc_genoprob(grav2, step=1, error_prob=0.002)

calc_genoprob <-
function(cross, step=0, off_end=0, stepwidth=c("fixed", "max"), pseudomarker_map=NULL,
         error_prob=1e-4, map_function=c("haldane", "kosambi", "c-f", "morgan"),
         lowmem=FALSE, quiet=TRUE, cores=1)
{
    # check inputs
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    if(error_prob < 0)
        stop("error_prob must be > 0")
    map_function <- match.arg(map_function)
    stepwidth <- match.arg(stepwidth)

    if(!lowmem) # use other version
        return(calc_genoprob2(cross=cross, step=step, off_end=off_end, stepwidth=stepwidth,
                              pseudomarker_map=pseudomarker_map, error_prob=error_prob,
                              map_function=map_function, quiet=quiet, cores=cores))

    # set up cluster; set quiet=TRUE if multi-core
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    # construct map at which to do the calculations
    # tolerance for matching marker and pseudomarker positions
    tol <- ifelse(step==0 || step>1, 0.01, step/100)
    # create the combined marker/pseudomarker map
    map <- insert_pseudomarkers(cross$gmap, step, off_end, stepwidth,
                                pseudomarker_map, tol)
    index <- map$index
    grid <- map$grid
    map <- map$map

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
    result <- list(probs=probs,
                   map=map,
                   crosstype=cross$crosstype,
                   is_x_chr=cross$is_x_chr,
                   is_female=cross$is_female,
                   cross_info=cross$cross_info,
                   alleles=cross$alleles,
                   alleleprobs=FALSE,
                   step=step,
                   off_end=off_end,
                   stepwidth=stepwidth,
                   error_prob=error_prob,
                   map_function=map_function)
    result$grid <- grid # include only if not NULL

    class(result) <- c("calc_genoprob", "list")

    result
}

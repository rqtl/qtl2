#' Guess phase of imputed genotypes
#'
#' Turn imputed genotypes into phased genotypes along chromosomes by
#' attempting to pick the phase that leads to the fewest recombination
#' events.
#'
#' @param cross Object of class \code{"cross2"}. For details, see the
#' \href{http://kbroman.org/qtl2/assets/vignettes/developer_guide.html}{R/qtl2 developer guide}.
#' @param geno Imputed genotypes, as a list of matrices, as from \code{\link{maxmarg}}.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If \code{0}, use \code{\link[parallel]{detectCores}}.)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by \code{\link[parallel]{makeCluster}}.
#'
#' @return
#'     If input cross is phase-known (e.g., recombinant inbred lines),
#'     the output will be the input \code{geno}. Otherwise, the output
#'     will be a list of three-dimensional arrays of imputed
#'     genotypes, individual x position x haplotype (1/2).
#'
#' @details We randonly assign the pair of alleles at the first locus
#'     to two haplotypes, and then work left to right, assigning
#'     alleles to haplotypes one locus at a time seeking the fewest
#'     recombination events. The results are subject to arbitrary and
#'     random choices. For example, to the right of a homozygous
#'     region, either orientation is equally reasonable.
#'
#' @export
#' @keywords utilities
#' @seealso \code{\link{maxmarg}}
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' gmap <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, gmap, error_prob=0.002)
#' imp_geno <- maxmarg(probs)
#' ph_geno <- guess_phase(iron, imp_geno)
guess_phase <-
    function(cross, geno, cores=1)
{
    if(n_chr(cross) != length(geno) ||
       names(cross$geno) != names(geno))
        stop("Mismatch in chromosomes between cross and geno")

    if(is_phase_known(cross)) return(geno)

    else {
        # set up cluster; use quiet=TRUE
        cores <- setup_cluster(cores, TRUE)

        is_x_chr <- cross$is_x_chr
        crosstype <- cross$crosstype
        is_female <- cross$is_female

        by_chr_func <- function(i) {
            if(crosstype=="f2") {
                if(is_x_chr[i]) result <- .guess_phase_f2X(t(geno[[i]]))
                else result <- .guess_phase_f2A(t(geno[[i]]))
            }
            else {
                if(is_x_chr[i]) result <- .guess_phase_X(t(geno[[i]]), crosstype, is_female)
                else result <- .guess_phase_A(t(geno[[i]]), crosstype)
            }
            dn <- dimnames(geno[[i]])
            result <- aperm(result, c(3,2,1))
            dimnames(result) <- list(dn[[1]], dn[[2]], c("mom", "dad"))
            return(result)
        }

        result <- cluster_lapply(cores, seq_along(geno), by_chr_func)
        names(result) <- names(geno)
        return(result)
    }

}

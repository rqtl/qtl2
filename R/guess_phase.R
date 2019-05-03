#' Guess phase of imputed genotypes
#'
#' Turn imputed genotypes into phased genotypes along chromosomes by
#' attempting to pick the phase that leads to the fewest recombination
#' events.
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param geno Imputed genotypes, as a list of matrices, as from [maxmarg()].
#' @param deterministic If TRUE, preferentially put smaller allele first when there's uncertainty.
#' If FALSE, the order of alleles is random in such cases.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return
#' If input cross is phase-known (e.g., recombinant inbred lines),
#' the output will be the input `geno`. Otherwise, the output
#' will be a list of three-dimensional arrays of imputed
#' genotypes, individual x position x haplotype (1/2).
#'
#' @details We randomly assign the pair of alleles at the first locus
#'     to two haplotypes, and then work left to right, assigning
#'     alleles to haplotypes one locus at a time seeking the fewest
#'     recombination events. The results are subject to arbitrary and
#'     random choices. For example, to the right of a homozygous
#'     region, either orientation is equally reasonable.
#'
#' @export
#' @keywords utilities
#' @seealso [maxmarg()]
#'
#' @examples
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' gmap <- insert_pseudomarkers(iron$gmap, step=1)
#' probs <- calc_genoprob(iron, gmap, error_prob=0.002)
#' imp_geno <- maxmarg(probs)
#' ph_geno <- guess_phase(iron, imp_geno)
guess_phase <-
    function(cross, geno, deterministic=FALSE, cores=1)
{
    # if cross is phase-known, don't change geno
    if(!is.cross2(cross))
        stop('Input cross must have class "cross2"')
    if(is_phase_known(cross)) return(geno)

    geno <- unclass(geno) # treat geno as a plain list
    if(!is.list(geno)) stop("geno should be a list of genotype matrices")

    # match chromosomes
    if(n_chr(cross) != length(geno) ||
       any(names(cross$geno) != names(geno))) {
        chr_cross <- chr_names(cross)
        chr_geno <- names(geno)
        if(!any(chr_cross %in% chr_geno))
            stop("No chromosomes in common between cross and geno")
        common_chr <- chr_cross[chr_cross %in% chr_geno]
        cross <- subset(cross, chr=common_chr)
        geno <- geno[common_chr]
    }

    # match individuals, if X chr
    if(any(cross$is_x_chr)) {
        ind_cross <- names(cross$is_female)
        xchr <- chr_names(cross)[cross$is_x_chr]
        ind_geno <- rownames(geno[[xchr[1]]])

        common_ind <- ind_cross[ind_cross %in% ind_geno]
        if(length(common_ind) == 0)
            stop("No individuals in common between cross and geno")
        cross <- subset(cross, ind=common_ind)
        for(i in names(geno)[cross$is_x_chr])
            geno[[i]] <- geno[[i]][common_ind,,drop=FALSE]
        if(length(common_ind) < length(ind_cross) ||
           length(common_ind) < length(ind_geno))
            warning("On X chr, considering only the ", length(common_ind),
                    ifelse(length(common_ind)==1, " individual", " individuals"),
                    " in common between cross and geno")
    }

    # set up cluster; use quiet=TRUE
    cores <- setup_cluster(cores, TRUE)

    is_x_chr <- cross$is_x_chr
    crosstype <- cross$crosstype
    is_female <- cross$is_female

    by_chr_func <- function(i) {
        if(crosstype=="f2") {
            if(is_x_chr[i]) result <- .guess_phase_f2X(t(geno[[i]]), deterministic)
            else result <- .guess_phase_f2A(t(geno[[i]]), deterministic)
        }
        else {
            if(is_x_chr[i]) result <- .guess_phase_X(t(geno[[i]]), crosstype, is_female, deterministic)
            else result <- .guess_phase_A(t(geno[[i]]), crosstype, deterministic)
        }
        dn <- dimnames(geno[[i]])
        result <- aperm(result, c(3,2,1))
        dimnames(result) <- list(dn[[1]], dn[[2]], c("mom", "dad"))
        return(result)
    }

    result <- cluster_lapply(cores, seq_along(geno), by_chr_func)
    names(result) <- names(geno)
    result
}

#' Find IBD segments for a set of strains
#'
#' Find IBD segments (regions with a lot of shared SNP genotypes) for
#' a set of strains
#'
#' @param geno List of matrices of founder genotypes. The matrices
#'     correspond to the genotypes on chromosomes and are arrayed as
#'     founders x markers.
#' @param map List of vectors of marker positions
#' @param min_lod Threshold for minimum LOD score for a segment
#' @param error_prob Genotyping error/mutation probability
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A data frame whose rows are IBD segments and whose columns
#'     are:
#' * Strain 1
#' * Strain 2
#' * Chromosome
#' * Left marker
#' * Right marker
#' * Left position
#' * Right position
#' * Left marker index
#' * Right marker index
#' * Interval length
#' * Number of markers
#' * Number of mismatches
#' * LOD score
#'
#' @details For each strain pair on each chromosome, we consider all
#' marker intervals and calculate a LOD score comparing the two
#' hypotheses: that the strains are IBD in the interval, vs. that
#' they are not. We assume that the two strains are homozygous at
#' all markers, and use the model from Broman and Weber (1999),
#' which assumes linkage equilibrium between markers and uses a
#' simple model for genotype frequencies in the presence of genotyping
#' errors or mutations.
#'
#' Note that inference of IBD segments is heavily dependent on how
#' SNPs were chosen to be genotyped. (For example, were the SNPs ascertained
#' based on their polymorphism between a particular strain pair?)
#'
#' @references
#' Broman KW, Weber JL (1999) Long homozygous chromosomal segments in
#' reference families from the Centre d’Étude du Polymorphisme Humain.
#' Am J Hum Genet 65:1493--1500.
#'
#' @examples
#' \dontrun{
#' # load DO data from Recla et al. (2014) Mamm Genome 25:211-222.
#' recla <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/master/DO_Recla/recla.zip")
#'
#' # grab founder genotypes and physical map
#' fg <- recla$founder_geno
#' pmap <- recla$pmap
#'
#' # find shared segments
#' (segs <- find_ibd_segments(fg, pmap, min_lod=10, error_prob=0.0001))
#' }
#'
#' @export
find_ibd_segments <-
    function(geno, map, min_lod=15, error_prob=0.001, cores=1)
{
    if(is.null(geno)) stop("geno is NULL")
    if(is.null(map)) stop("map is NULL")

    if(!is_number(min_lod)) stop("min_lod should be a single number")

    n_str <- vapply(geno, nrow, 1)
    if(length(unique(n_str)) != 1)
        stop("The geno matrices have different numbers of rows; they should be strains x markers")
    n_str <- n_str[1]
    if(n_str < 2)
        stop("Need at least two strains")

    nmar <- vapply(geno, ncol, 1)
    nmar_map <- vapply(map, length, 1)
    if(length(nmar) != length(nmar_map))
        stop("length(geno) != length(map)")
    if(!all(nmar == nmar_map))
        stop("geno and map have different numbers of markers")

    if(!is_pos_number(error_prob) || error_prob >= 1)
        stop("error_prob should be a single number in (0,1)")

    # get strain pairs
    str1 <- matrix(rep(1:n_str, n_str), ncol=n_str)
    str2 <- matrix(rep(1:n_str, n_str), ncol=n_str, byrow=TRUE)
    str_list <- data.frame(strain1=str1[upper.tri(str1)],
                           strain2=str2[upper.tri(str2)])
    str_list <- str_list[order(str_list[,1], str_list[,2]),]
    n_pair <- nrow(str_list)

    n_chr <- length(geno)
    run_list <- data.frame(str1 = rep(str_list[,1], n_chr),
                           str2 = rep(str_list[,2], n_chr),
                           chr = rep(1:n_chr, each=n_pair))

    freq <- lapply(geno, function(a) colMeans(a==1, na.rm=TRUE))

    str_names <- rownames(geno[[1]])
    if(is.null(str_names)) {
        if(n_str <= 26) str_names <- LETTERS[1:n_str]
        else str_names <- as.character(1:n_str)
    }

    by_run_func <- function(i) {
        chr <- run_list$chr[i]
        chrnam <- names(map)[chr]
        str1 <- run_list$str1[i]
        str2 <- run_list$str2[i]

        g1 <- geno[[chr]][str1,]
        g2 <- geno[[chr]][str2,]
        m <- map[[chr]]
        p <- freq[[chr]]

        keep <- which(!is.na(g1) & !is.na(g2) & p > 0 & p < 1)

        result <- .find_ibd_segments(g1[keep], g2[keep], p[keep], error_prob)
        result <- as.data.frame(result)

        # drop rows (low LOD or overlapping interval with higher LOD)
        result <- result[result[,3] >= min_lod & result[,6] > 0.5,,drop=FALSE]

        if(nrow(result)==0) return(NULL)

        data.frame(strain1=str_names[str1],
                   strain2=str_names[str2],
                   chr=chrnam,
                   left_marker=names(m)[result[,1]],
                   right_marker=names(m)[result[,2]],
                   left_pos=m[result[,1]],
                   right_pos=m[result[,2]],
                   left_index=as.integer(keep[result[,1]]),
                   right_index=as.integer(keep[result[,2]]),
                   int_length=m[result[,2]] - m[result[,1]],
                   n_mar = as.integer(result[,2] - result[,1] + 1),
                   n_mismatch=as.integer(result[,5]),
                   lod=result[,3],
                   stringsAsFactors=FALSE)

    }

    # set up cluster; set quiet=TRUE if multi-core
    quiet <- TRUE
    cores <- setup_cluster(cores, quiet)
    if(!quiet && n_cores(cores)>1) {
        message(" - Using ", n_cores(cores), " cores")
        quiet <- TRUE # make the rest quiet
    }

    result_list <- cluster_lapply(cores, 1:nrow(run_list), by_run_func)
    if(all(vapply(result_list, is.null, TRUE))) return(NULL) # no segments at all!
    result_list <- result_list[!vapply(result_list, is.null, TRUE)] # just save runs that returned some rows

    # combine the results into a single data frame
    result <- do.call("rbind", result_list)

    # reorder by chromosome, left pos, right pos, strain 1, strain 2
    result <- result[order(factor(result$chr, names(map)),
                           result$left_pos,
                           result$right_pos,
                           factor(result$strain1, str_names),
                           factor(result$strain2, str_names)), ]

    rownames(result) <- as.character(1:nrow(result))
    result
}

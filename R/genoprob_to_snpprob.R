#' Convert genotype probabilities to SNP probabilities
#'
#' For multi-parent populations, convert use founder genotypes at a
#' set of SNPs to convert founder-based genotype probabilities to SNP
#' genotype probabilities.
#'
#' @param genoprobs Genotype probabilities as
#' calculated by \code{\link[qtl2geno]{calc_genoprob}}.
#' @param map Physical map of markers and pseudomarkers in the
#' \code{genoprobs} object.
#' @param snpinfo Data frame with SNP information with the following columns:
#' \itemize{
#' \item \code{chr} - Character string or factor with chromosome
#' \item \code{pos} - Position (in same units as in the \code{"map"}
#'     attribute in \code{genoprobs}.
#' \item \code{sdp} - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' \item \code{snp} - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' }
#' @param tol Tolerance for determining whether a SNP is exactly at a
#' position at which genotype probabilities were already calculated.
#'
#' @return A list with two components:
#' \itemize{
#' \item \code{probs} - An object like the \code{genoprobs} input,
#' but with imputed genotype probabilities at the selected SNPs.
#' An added attribute \code{"snpinfo"} contains the input \code{snpinfo} with
#' an additional column, \code{index}, that is an index to the
#' genotype probabilities. (In cases where multiple SNPs are contained
#' in the same interval and have the same strain distribution pattern
#' (SDP), the probabilities would be all the same and so the results
#' contain the probabilities for just one of them.)
#' \item \code{map} - An object like the \code{map} input, but
#' including the positions of imputed SNPs.
#' }
#'
#' If the input \code{genoprobs} is for allele probabilities, the
#' \code{probs} output has just two probability columns (for the two SNP
#' alleles). If the input has a full set of \eqn{n(n+1)/2}
#' probabilities for \eqn{n} strains, the \code{probs} output has 3 probabilities
#' (for the three SNP genotypes). If the input has full genotype
#' probabilities for the X chromosome (\eqn{n(n+1)/2} genotypes for
#' the females followed by \eqn{n} hemizygous genotypes for the
#' males), the output has 5 probabilities: the 3 female SNP genotypes
#' followed by the two male hemizygous SNP genotypes.
#'
#' @details We first split the SNPs by chromosome and identify the
#' intervals in the genotype probabilities that contain each. For SNPs
#' within \code{tol} of a position at which the genotype probabilities
#' were calculated, we use those probabilities. For SNPs contained
#' within an interval, we use the average of the probabilities for the
#' two endpoints. We then collapse the probabilities according to the
#' strain distribution pattern.
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' library(qtl2geno)
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DO_Recla/recla.zip")
#' recla <- read_cross2(file)
#' recla <- recla[c(1:2,53:54), c("19","X")] # subset to 4 mice and 2 chromosomes
#' probs <- calc_genoprob(recla, error_prob=0.002)
#'
#' # founder genotypes for a set of SNPs
#' snpgeno <- rbind(m1=c(3,1,1,3,1,1,1,1),
#'                  m2=c(1,3,1,3,1,3,1,3),
#'                  m3=c(1,1,1,1,3,3,3,3),
#'                  m4=c(1,3,1,3,1,3,1,3))
#' sdp <- calc_sdp(snpgeno)
#' snpinfo <- data.frame(chr=c("19", "19", "X", "X"),
#'                       pos=c(40.36, 40.53, 110.91, 111.21),
#'                       sdp=sdp,
#'                       snp=c("m1", "m2", "m3", "m4"), stringsAsFactors=FALSE)
#'
#' # collapse to SNP genotype probabilities
#' snpprobs <- genoprob_to_snpprob(probs, recla$pmap, snpinfo)
#'
#' # could also first convert to allele probs
#' aprobs <- genoprob_to_alleleprob(probs)
#' snpaprobs <- genoprob_to_snpprob(aprobs, recla$pmap, snpinfo)
#' }
#'
#' @export
genoprob_to_snpprob <-
    function(genoprobs, map, snpinfo, tol=1e-8)
{
    uchr <- unique(snpinfo$chr)
    chrID <- names(genoprobs)
    if(!all(uchr %in% chrID)) {
        mischr <- uchr[!(uchr %in% chrID)]
        stop("Not all chr found in genoprobs: ", paste(mischr, collapse=","))
    }
    # reorder
    uchr <- factor(factor(uchr, levels=chrID))

    # if more than one chromosome:
    if(length(uchr) > 1) {
        uchr <- as.character(uchr) # convert back to character

        ### split snpinfo by chr, ordered as in genoprobs
        chr <- factor(snpinfo$chr, levels=uchr)
        snpinfo_spl <- split(snpinfo, chr)

        ### loop over chromosomes, using recursion
        empty_vec <- vector("list", length(uchr))
        names(empty_vec) <- uchr
        results <- list(probs=empty_vec, map=empty_vec)
        attr(results$probs, "snpinfo") <- empty_vec
        for(i in seq(along=uchr)) {
            tmp <- genoprob_to_snpprob(genoprobs, map, snpinfo_spl[[i]])
            results$probs[[i]] <- tmp$probs[[1]] # just the array
            attr(results$probs, "snpinfo")[[i]] <- attr(tmp$probs, "snpinfo")[[1]]
            results$map[[i]] <- tmp$map[[1]]
        }

        ### add attributes
        attr(results$probs, "crosstype") <- "snps"
        attr(results$probs, "alleles") <- c("A", "B")
        attr(results$map, "is_x_chr") <- attr(results$probs, "is_x_chr") <-
            attr(genoprobs, "is_x_chr")[uchr]
        class(results$probs) <- c("calc_genoprob", "list")

        return(results)
    }

    # one chromosome

    # make chromosome a character string again
    uchr <- as.character(uchr)

    ### number of alleles
    alleles <- attr(genoprobs, "alleles")
    n_alleles <- length(alleles)

    ### find snps in map
    this_map <- map[[uchr]]
    snploc <- find_intervals(snpinfo$pos, this_map, tol)
    interval <- snploc[,1]
    on_map <- (snploc[,2]==1)

    # drop snps outside of range
    snps2drop <- (interval < 0 | (interval >= length(this_map)-1 & !on_map) |
                  snpinfo$sdp<1 | snpinfo$sdp >(2^n_alleles-2))
    if(any(snps2drop)) {
        snpinfo <- snpinfo[!snps2drop,,drop=FALSE]
        interval <- interval[!snps2drop]
        on_map <- interval[!snps2drop]
    }
    if(nrow(snpinfo) == 0)
        stop("No SNPs within range")

    ### find unique (interval, on_map, sdp) patterns
    pat <- paste(interval, on_map, snpinfo$sdp, sep=":")
    upat <- unique(pat)
    ### order unique patterns by position
    pos <- snpinfo$pos[match(upat, pat)]
    upat <- upat[order(pos)]
    snpinfo$index <- match(pat, upat)
    ### first row with each unique pattern
    urow <- (1:nrow(snpinfo))[match(upat, pat)]
    ### distinct bits
    interval <- interval[urow]
    on_map <- on_map[urow]
    sdp <- snpinfo$sdp[urow]

    # subset to the single chromosome
    empty_vec <- vector("list", 1)
    names(empty_vec) <- uchr
    results <- list(probs=empty_vec, map=empty_vec)
    attr(results$probs, "snpinfo") <- empty_vec

    ### genoprobs -> SNP genotype probs
    if(ncol(genoprobs[[uchr]]) == n_alleles) { # allele probs
        results$probs[[1]] <- .alleleprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("A", "B")
    }
    else if(ncol(genoprobs[[uchr]]) == n_alleles*(n_alleles+1)/2) { # autosomal genotype probs
        results$probs[[1]] <- .genoprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("AA", "AB", "BB")
    }
    else if(ncol(genoprobs[[uchr]]) == n_alleles + n_alleles*(n_alleles+1)/2) { # X chr genotype probs
        results$probs[[1]] <- .Xgenoprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("AA", "AB", "BB", "AY", "BY")
    }
    else {
        stop("genoprobs has ", ncol(genoprobs[[uchr]]),
             " columns but there are ", n_alleles, " alleles")
    }

    # add dimnames
    snpnames <- snpinfo$snp[urow]
    if(is.null(snpnames)) snpnames <- rownames(snpinfo)[urow]
    dimnames(results$probs[[1]]) <- list(rownames(genoprobs[[uchr]]),
                                        coln, snpnames)

    ### add attributes
    tmp <- list("1"=snpinfo)
    names(tmp) <- uchr
    attr(results$probs, "snpinfo") <- tmp

    snpmap <- list(snpinfo$pos[urow])
    names(snpmap[[1]]) <- snpnames
    names(snpmap) <- uchr
    results$map <- snpmap
    attr(results$map, "is_x_chr") <- attr(results$probs, "is_x_chr") <-
        attr(genoprobs, "is_x_chr")[uchr]
    attr(results$probs, "crosstype") <- "snps"
    attr(results$probs, "alleles") <- c("A", "B")
    class(results$probs) <- c("calc_genoprob", "list")

    snpinfo <- list("1"=snpinfo)
    names(snpinfo) <- uchr
    attr(results$probs, "snpinfo") <- snpinfo

    results
}

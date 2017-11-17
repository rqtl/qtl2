#' Convert genotype probabilities to SNP probabilities
#'
#' For multi-parent populations, convert use founder genotypes at a
#' set of SNPs to convert founder-based genotype probabilities to SNP
#' genotype probabilities.
#'
#' @md
#'
#' @param genoprobs Genotype probabilities as
#' calculated by [qtl2geno::calc_genoprob()].
#'
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived with
#'     [index_snps()]):
#' * `chr` - Character string or factor with chromosome
#' * `pos` - Position (in same units as in the `"map"`
#'     attribute in `genoprobs`.
#' * `sdp` - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' * `snp` - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' * `index` - Indices that indicate equivalent
#'     groups of SNPs, calculated by [index_snps()].
#' * `intervals` - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' * `on_map` - Indicate whether SNP coincides with a marker
#'     in the `genoprobs`
#'
#' @return An object like the `genoprobs` input, but with imputed
#' genotype probabilities at the selected SNPs indicated in
#' `snpinfo$index`.
#'
#' If the input `genoprobs` is for allele probabilities, the
#' `probs` output has just two probability columns (for the two SNP
#' alleles). If the input has a full set of \eqn{n(n+1)/2}
#' probabilities for \eqn{n} strains, the `probs` output has 3 probabilities
#' (for the three SNP genotypes). If the input has full genotype
#' probabilities for the X chromosome (\eqn{n(n+1)/2} genotypes for
#' the females followed by \eqn{n} hemizygous genotypes for the
#' males), the output has 5 probabilities: the 3 female SNP genotypes
#' followed by the two male hemizygous SNP genotypes.
#'
#' @details We first split the SNPs by chromosome and use
#' `snpinfo$index` to subset to non-equivalent SNPs.
#' `snpinfo$interval` indicates the intervals in the genotype
#' probabilities that contain each. For SNPs contained within an
#' interval, we use the average of the probabilities for the two
#' endpoints. We then collapse the probabilities according to the
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
#' # identify groups of equivalent SNPs
#' snpinfo <- index_snps(recla$pmap, snpinfo)
#'
#' # collapse to SNP genotype probabilities
#' snpprobs <- genoprob_to_snpprob(probs, snpinfo)
#'
#' # could also first convert to allele probs
#' aprobs <- genoprob_to_alleleprob(probs)
#' snpaprobs <- genoprob_to_snpprob(aprobs, snpinfo)
#' }
#'
#' @seealso [index_snps()], [qtl2geno::calc_genoprob()]
#' @export
genoprob_to_snpprob <-
    function(genoprobs, snpinfo)
{
    if(nrow(snpinfo)==0) {
        result <- genoprobs[,names(genoprobs)[1]]
        if(length(attr(result, "alleles")) == ncol(result[[1]])) {
            result[[1]] <- result[[1]][,1:2,numeric(0),drop=FALSE]
            colnames(result[[1]]) <- c("A", "B")
        } else if(attr(result, "is_x_chr")[1]) {
            result[[1]] <- result[[1]][,1:5,numeric(0),drop=FALSE]
            colnames(result[[1]]) <- c("AA", "AB", "BB", "AY", "BY")
        } else {
            result[[1]] <- result[[1]][,1:3,numeric(0),drop=FALSE]
            colnames(result[[1]]) <- c("AA", "AB", "BB")
        }
        return(result)
    }


    uchr <- unique(snpinfo$chr)
    chrID <- names(genoprobs)
    if(!all(uchr %in% chrID)) {
        mischr <- uchr[!(uchr %in% chrID)]
        stop("Not all chr found in genoprobs: ", paste(mischr, collapse=","))
    }
    # reorder
    uchr <- factor(factor(uchr, levels=chrID))

    # check for index, interval, and on_map
    if(!all(c("index", "interval", "on_map") %in% colnames(snpinfo)))
        stop('snpinfo should contain columns "index", "interval", and "on_map".')

    # if more than one chromosome:
    if(length(uchr) > 1) {
        uchr <- as.character(uchr) # convert back to character

        ### split snpinfo by chr, ordered as in genoprobs
        chr <- factor(snpinfo$chr, levels=uchr)
        snpinfo_spl <- split(snpinfo, chr)

        ### loop over chromosomes, using recursion
        results <- vector("list", length(uchr))
        names(results) <- uchr
        for(i in seq_along(uchr))
            results[[i]] <- genoprob_to_snpprob(genoprobs, snpinfo_spl[[i]])[[1]]

        ### add attributes
        attr(results, "crosstype") <- "snps"
        attr(results, "alleles") <- c("A", "B")
        attr(results, "is_x_chr") <- attr(genoprobs, "is_x_chr")[uchr]
        class(results) <- c("calc_genoprob", "list")

        return(results)
    }

    # one chromosome

    # make chromosome a character string again
    uchr <- as.character(uchr)

    ### number of alleles
    alleles <- attr(genoprobs, "alleles")
    n_alleles <- length(alleles)

    ### subset the snpinfo
    snpinfo <- snpinfo[sort(unique(snpinfo$index)),,drop=FALSE]

    ### grab stuff
    sdp <- snpinfo$sdp
    interval <- snpinfo$interval
    on_map <- snpinfo$on_map

    if(any(interval < 0 | interval > dim(genoprobs[[uchr]])[3]-1))
        stop("interval values outside of the range [0, ",
             dim(genoprobs[[uchr]])[3]-1, "].")

    if(any(sdp < 1 | sdp > 2^n_alleles-1))
        stop("sdp values outside of the range [1, ",
             2^n_alleles-1, "].")

    # subset to the single chromosome
    results <- vector("list", 1)
    names(results) <- uchr

    ### genoprobs -> SNP genotype probs
    if(ncol(genoprobs[[uchr]]) == n_alleles) { # allele probs
        results[[1]] <- .alleleprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("A", "B")
    }
    else if(ncol(genoprobs[[uchr]]) == n_alleles*(n_alleles+1)/2) { # autosomal genotype probs
        results[[1]] <- .genoprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("AA", "AB", "BB")
    }
    else if(ncol(genoprobs[[uchr]]) == n_alleles + n_alleles*(n_alleles+1)/2) { # X chr genotype probs
        results[[1]] <- .Xgenoprob_to_snpprob(genoprobs[[uchr]], sdp, interval, on_map)
        coln <- c("AA", "AB", "BB", "AY", "BY")
    }
    else {
        stop("genoprobs has ", ncol(genoprobs[[uchr]]),
             " columns but there are ", n_alleles, " alleles")
    }

    # add dimnames
    snpnames <- snpinfo$snp
    if(is.null(snpnames)) snpnames <- rownames(snpinfo)
    dimnames(results[[1]]) <- list(rownames(genoprobs[[uchr]]),
                                        coln, snpnames)

    attr(results, "is_x_chr") <- attr(genoprobs, "is_x_chr")[uchr]
    attr(results, "crosstype") <- "snps"
    attr(results, "alleles") <- c("A", "B")
    class(results) <- c("calc_genoprob", "list")

    results
}

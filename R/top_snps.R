#' Create table of top snp associations
#'
#' Create a table of the top snp associations
#'
#' @param scan1_output Output of \code{\link[qtl2scan]{scan1}}.
#' Should contain a component \code{"snpinfo"}, as when
#' \code{\link[qtl2scan]{scan1}} is run with SNP probabilities
#' produced by \code{\link[qtl2scan]{genoprob_to_snpprob}}.
#'
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived from with
#'     \code{\link{index_snps}}):
#' \itemize{
#' \item \code{chr} - Character string or factor with chromosome
#' \item \code{pos} - Position (in same units as in the \code{"map"}
#'     attribute in \code{genoprobs}.
#' \item \code{sdp} - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' \item \code{snp} - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' \item \code{index} - Indices that indicate equivalent
#'     groups of SNPs, calculated by \code{\link{index_snps}}.
#' \item \code{intervals} - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' \item \code{on_map} - Indicate whether SNP coincides with a marker
#'     in the \code{genoprobs}
#' }
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param drop Show all SNPs with LOD score within this amount of the
#' maximum SNP association.
#'
#' @examples
#' \dontrun{
#' # load example DO data from web
#' library(qtl2geno)
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#'
#' # subset to chr 2
#' DOex <- DOex[,"2"]
#'
#' # calculate genotype probabilities and convert to allele probabilities
#' pr <- calc_genoprob(DOex, error_prob=0.002)
#' apr <- genoprob_to_alleleprob(pr)
#'
#' # download snp info from web
#' tmpfile <- tempfile()
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/c2_snpinfo.rds")
#' download.file(file, tmpfile, quiet=TRUE)
#' snpinfo <- readRDS(tmpfile)
#' unlink(tmpfile)
#'
#' # calculate strain distribution patterns
#' snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
#'
#' # identify groups of equivalent SNPs
#' snpinfo <- index_snps(DOex$pmap, snpinfo)
#'
#' # convert to snp probabilities
#' snppr <- genoprob_to_snpprob(apr, snpinfo)
#'
#' # perform SNP association analysis (here, ignoring residual kinship)
#' out_snps <- scan1(snppr, DOex$pheno)
#'
#' # table with top SNPs
#' top_snps(out_snps, snpinfo)
#'
#' # top SNPs among the distinct subset at which calculations were performed
#' top_snps(out_snps, snpinfo, show_all_snps=FALSE)
#'
#' # top SNPs within 0.5 LOD of max
#' top_snps(out_snps, snpinfo, 0.5)
#' }
#' @export
#' @seealso \code{\link{index_snps}}, \code{\link{genoprob_to_snpprob}}, \code{\link[qtl2plot]{plot_snpasso}}
top_snps <-
    function(scan1_output, snpinfo, drop=1.5, show_all_snps=TRUE)
{
    uindex <- unique(snpinfo$index)
    if(length(uindex) != nrow(scan1_output))
        stop("Something is wrong with snpinfo$index.\n",
             "      length(unique(snpinfo$index)) [",
             length(unique(snpinfo$index)), "] != nrow(scan1_output) [",
             length(lod), "].")

    if(any(snpinfo$index[uindex] != uindex))
        stop("Something is wrong with snpinfo$index.\n",
             "      snpinfo$index[u] should == u for values in snpinfo$index")

    map <- snpinfo_to_map(snpinfo)

    chr <- names(map)
    if(length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    # deal with possibly > 1 chr
    lod <- unclass(scan1_output[map, chr, ])

    map <- map[[chr]]
    snpinfo <- snpinfo[snpinfo$chr==chr,,drop=FALSE]

    if(ncol(scan1_output) > 1)
        warning("Considering only the first LOD score column")

    keep <- which(!is.na(lod) & lod >= max(lod, na.rm=TRUE) - drop)

    # reverse the snp index
    revindex <- rev_snp_index(snpinfo)

   if(show_all_snps) { # expand to all related SNPs
        snpinfo <- snpinfo[revindex %in% keep,,drop=FALSE]
        revindex <- revindex[revindex %in% keep]
        snpinfo$lod <- lod[revindex]
    } else { # just keep the SNPs that were used
        snpinfo$lod <- lod[revindex]
        snpinfo <- snpinfo[snpinfo$snp %in% rownames(scan1_output)[keep],]
    }

    snpinfo
}

# snpinfo to map
snpinfo_to_map <-
    function(snpinfo)
{
    uindex <- sort(unique(snpinfo$index))
    if(any(snpinfo$index < 1 | snpinfo$index > nrow(snpinfo)))
        stop("snpinfo$index values outside of range [1, ",
             nrow(snpinfo), "]")

    uchr <- unique(snpinfo$chr)
    chr <- factor(snpinfo$chr, levels=uchr)

    map <- split(snpinfo$pos, chr)
    snp <- split(snpinfo$snp, chr)
    index <- split(snpinfo$index, chr)
    for(i in seq(along=map)) {
        u <- unique(index[[i]])
        map[[i]] <- map[[i]][u]
        names(map[[i]]) <- snp[[i]][u]
    }

    names(map) <- uchr

    map
}

# reverse index
rev_snp_index <-
    function(snpinfo)
{
    index_spl <- split(1:nrow(snpinfo), snpinfo$index)
    revindex <- rep(seq(along=index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex

    revindex
}

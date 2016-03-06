#' Create table of top snp associations
#'
#' Create a table of the top snp associations
#'
#' @param scan1output Output of \code{\link[qtl2scan]{scan1}} or
#' \code{\link[qtl2scan]{scan1_lmm}}. Should contain an attribute,
#' \code{"snpinfo"}, as when \code{\link[qtl2scan]{scan1}} or
#' \code{\link[qtl2scan]{scan1_lmm}} are run with SNP probabilities
#' produced by \code{\link[qtl2scan]{genoprob_to_snpprob}}.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param drop Show all SNPs with LOD score within this amount of the
#' maximum SNP association.
#'
#' @export
#' @seealso \code{\link{genoprob_to_snpprob}}, \code{\link[qtl2plot]{plot_snpasso}}
top_snps <-
    function(scan1output, drop=1.5, show_all_snps=TRUE)
{
    map <- attr(scan1output, "map")
    snpinfo <- attr(scan1output, "snpinfo")
    if(is.null(snpinfo)) stop("No snpinfo found")
    if(is.null(map)) stop("No map found")

    chr <- names(map)
    if(length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }
    map <- map[[chr]]
    snpinfo <- snpinfo[[chr]]

    if(ncol(scan1output) > 1)
        warning("Considering only the first LOD score column")

    lod <- scan1output[seq(along=map),1]
    keep <- which(!is.na(lod) & lod > max(lod, na.rm=TRUE) - drop)

    if(show_all_snps) { # expand to all related SNPs
        snpinfo <- snpinfo[snpinfo$index %in% keep,]
        snpinfo$lod <- lod[snpinfo$index]
    } else { # just keep the SNPs that were used
        snpinfo$lod <- lod[snpinfo$index]
        snpinfo <- snpinfo[snpinfo$snp %in% rownames(scan1output)[keep],]
    }

    snpinfo
}

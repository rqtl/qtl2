#' Create table of top snp associations
#'
#' Create a table of the top snp associations
#'
#' @param scan1_output Output of \code{\link[qtl2scan]{scan1}}.
#' Should contain a component \code{"snpinfo"}, as when
#' \code{\link[qtl2scan]{scan1}} is run with SNP probabilities
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
    function(scan1_output, drop=1.5, show_all_snps=TRUE)
{
    map <- scan1_output$map
    if(is.null(map)) stop("No map found")
    snpinfo <- scan1_output$snpinfo
    if(is.null(snpinfo)) stop("No snpinfo found")

    chr <- names(map)
    if(length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    # deal with possibly > 1 chr
    lod <- subset(scan1_output, chr=chr)$lod

    map <- map[[chr]]
    snpinfo <- snpinfo[[chr]]

    if(ncol(scan1_output$lod) > 1)
        warning("Considering only the first LOD score column")

    keep <- which(!is.na(lod) & lod > max(lod, na.rm=TRUE) - drop)

    if(show_all_snps) { # expand to all related SNPs
        snpinfo <- snpinfo[snpinfo$index %in% keep,]
        snpinfo$lod <- lod[snpinfo$index]
    } else { # just keep the SNPs that were used
        snpinfo$lod <- lod[snpinfo$index]
        snpinfo <- snpinfo[snpinfo$snp %in% rownames(scan1_output$lod)[keep],]
    }

    snpinfo
}

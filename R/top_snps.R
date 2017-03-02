#' Create table of top snp associations
#'
#' Create a table of the top snp associations
#'
#' @param scan1_output Output of \code{\link[qtl2scan]{scan1}}.
#' Should contain a component \code{"snpinfo"}, as when
#' \code{\link[qtl2scan]{scan1}} is run with SNP probabilities
#' produced by \code{\link[qtl2scan]{genoprob_to_snpprob}}.

#' @param map A list of vectors of SNP locations; generally this would
#' be the \code{"map"} component of the result produced by
#' \code{\link{genoprob_to_snpprob}}.
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
#' # convert to snp probabilities
#' snppr <- genoprob_to_snpprob(apr, DOex$pmap, snpinfo)
#'
#' # perform SNP association analysis (here, ignoring residual kinship)
#' out_snps <- scan1(snppr$probs, DOex$pheno)
#'
#' # table with top SNPs
#' top_snps(out_snps, snppr$map)
#'
#' # top SNPs among the distinct subset at which calculations were performed
#' top_snps(out_snps, snppr$map, show_all_snps=FALSE)
#'
#' # top SNPs within 0.5 LOD of max
#' top_snps(out_snps, snppr$map, 0.5)
#' }
#' @export
#' @seealso \code{\link{genoprob_to_snpprob}}, \code{\link[qtl2plot]{plot_snpasso}}
top_snps <-
    function(scan1_output, map, drop=1.5, show_all_snps=TRUE)
{
    if(is.null(map)) stop("No map found")
    snpinfo <- attr(scan1_output, "snpinfo")
    if(is.null(snpinfo)) stop("No snpinfo found")

    chr <- names(map)
    if(length(chr) > 1) {
        warning("Considering only chromosome ", chr)
        chr <- chr[1]
    }

    # deal with possibly > 1 chr
    lod <- unclass(scan1_output[map, chr, ])

    map <- map[[chr]]
    snpinfo <- snpinfo[[chr]]

    if(ncol(scan1_output) > 1)
        warning("Considering only the first LOD score column")

    keep <- which(!is.na(lod) & lod >= max(lod, na.rm=TRUE) - drop)

    if(show_all_snps) { # expand to all related SNPs
        snpinfo <- snpinfo[snpinfo$index %in% keep,]
        snpinfo$lod <- lod[snpinfo$index]
    } else { # just keep the SNPs that were used
        snpinfo$lod <- lod[snpinfo$index]
        snpinfo <- snpinfo[snpinfo$snp %in% rownames(scan1_output)[keep],]
    }

    snpinfo
}

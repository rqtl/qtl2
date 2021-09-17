#' Create snp information table for a cross
#'
#' Create a table of snp information from a cross, for use with [scan1snps()].
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#'
#' @return A data frame of SNP information with the following columns:
#' * `chr` - Character string or factor with chromosome
#' * `pos` - Position (in same units as in the `"map"`
#'     attribute in `genoprobs`.
#' * `snp` - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' * `sdp` - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' SNPs with missing founder genotypes are omitted.
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DO_Recla/recla.zip")
#' recla <- read_cross2(file)
#' snpinfo <- create_snpinfo(recla)
#'
#' # calculate genotype probabilities
#' pr <- calc_genoprob(recla, error_prob=0.002, map_function="c-f")
#'
#' # index the snp information
#' snpinfo <- index_snps(recla$pmap, snpinfo)
#'
#' # sex covariate
#' sex <- setNames((recla$covar$Sex=="female")*1, rownames(recla$covar))
#'
#' # perform a SNP scan
#' out <- scan1snps(pr, recla$pmap, recla$pheno[,"bw"], addcovar=sex, snpinfo=snpinfo)
#'
#' # plot the LOD scores
#' plot(out$lod, snpinfo, altcol="green3")
#' }
#'
#' @seealso [index_snps()], [scan1snps()], [genoprob_to_snpprob()]
#' @export
create_snpinfo <-
    function(cross)
{
    if(!is.cross2(cross)) {
        stop("Input should be a cross2 object")
    }

    if("pmap" %in% names(cross)) map <- cross$pmap
    else if("gmap" %in% names(cross)) map <- cross$gmap
    else stop("cross contains neither pmap nor gmap")

    # convert map to a data frame
    nmar <- vapply(map, length, 1) # no. markers per chromosome
    markers <- unlist(lapply(map, names))
    map <- data.frame(chr=rep(names(map), nmar),
                      pos=unlist(map),
                      snp=markers,
                      stringsAsFactors=FALSE)
    rownames(map) <- map$snp

    # founder genotypes -> SDP
    if(!("founder_geno" %in% names(cross))) {
        stop("cross does not contain founder_geno")
    }
    fg <- do.call("cbind", cross$founder_geno)

    # drop markers with missing founders or only one allele present in founders
    mar2drop <- (colSums(is.na(fg) | fg==0) > 0) |
        apply(fg, 2, function(g) length(unique(g[!is.na(g) & g!=0]))==1)

    # drop markers with missing data; calculate founder strain distribution patterns (SDPs)
    cbind(map[!mar2drop,,drop=FALSE],
          sdp=calc_sdp(t(fg[,!mar2drop])))
}

#' Recode SNPs by major allele
#'
#' For multi-parent populations with founder genotypes, recode the raw
#' SNP genotypes so that `1` means homozygous for the major allele in the
#' founders.
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#'
#' @return The input cross object with the raw SNP genotypes recoded so that
#' `1` is homozygous for the major alleles in the founders.
#'
#' @export
#' @keywords utilities
#' @seealso [calc_raw_founder_maf()], [calc_raw_maf()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' DOex <- recode_snps(DOex)
#' }
recode_snps <-
    function(cross)
{
    if(!("founder_geno" %in% names(cross))) {
        stop("cross does not contain founder genotypes")
    }

    for(chr in seq_along(cross$founder_geno)) {
        fg <- cross$founder_geno[[chr]]
        fg[fg!=1 & fg!=3] <- NA
        p1 <- colMeans(fg==1, na.rm=TRUE)
        recode <- (!is.na(p1) & p1 < 0.5)

        g <- cross$geno[[chr]]
        g[g==0] <- NA

        g[,recode] <- 4 - g[,recode]
        fg[,recode] <- 4 - fg[,recode]

        g[is.na(g)] <- 0
        fg[is.na(fg)] <- 0

        cross$founder_geno[[chr]] <- fg
        cross$geno[[chr]] <- g
    }

    cross
}

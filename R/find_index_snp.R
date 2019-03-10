#' Find name of indexed snp
#'
#' For a particular SNP, find the name of the corresponding indexed SNP.
#'
#' @param snpinfo Data frame with SNP information with the following columns:
#' * `chr` - Character string or factor with chromosome
#' * `index` - Numeric index of equivalent, indexed SNP, as produced by [index_snps()].
#' * `snp` - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' @param snp Name of snp to look for (can be a vector).
#'
#' @return A vector of SNP IDs (the corresponding indexed SNPs), with NA if a SNP is not found.
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DO_Recla/recla.zip")
#' recla <- read_cross2(file)
#'
#' # founder genotypes for a set of SNPs
#' snpgeno <- rbind(m1=c(3,1,1,3,1,1,1,1),
#'                  m2=c(3,1,1,3,1,1,1,1),
#'                  m3=c(1,1,1,1,3,3,3,3),
#'                  m4=c(1,3,1,3,1,3,1,3))
#' sdp <- calc_sdp(snpgeno)
#' snpinfo <- data.frame(chr=c("19", "19", "X", "X"),
#'                       pos=c(40.36, 40.53, 110.91, 111.21),
#'                       sdp=sdp,
#'                       snp=c("m1", "m2", "m3", "m4"), stringsAsFactors=FALSE)
#'
#' # update snp info by adding the SNP index column
#' snpinfo <- index_snps(recla$pmap, snpinfo)
#'
#' # find indexed snp for a particular snp
#' find_index_snp(snpinfo, "m3")
#' }
#'
#' @seealso [find_marker()]
#' @export
find_index_snp <-
    function(snpinfo, snp)
{
    if(!("index" %in% colnames(snpinfo)))
        stop('snpinfo does not contain an "index" column; run index_snps().')
    if(!("chr" %in% colnames(snpinfo)))
        stop('snpinfo does not contain an "chr" column.')

    if(length(unique(snpinfo$chr)) > 1) { # more than one chromosome
        # split the snpinfo table by chromosome
        snpinfo <- split(snpinfo, snpinfo$chr)

        # vector of NAs to contain the results
        result <- as.character(rep(NA, length(snp)))

        for(i in seq(along=snpinfo)) { # loop over chromosomes
            # if found everything already, just return...
            if(!any(is.na(result))) return(result)

            # grab snp names for this chromosome
            snps <- snpinfo[[i]]$snp
            if(is.null(snps)) snps <- rownames(snpinfo[[i]])
            index <- snpinfo[[i]]$index

            # indexed snps for those that have not yet been found
            result[is.na(result)] <- snps[index[match(snp[is.na(result)], snps)]]
        }

        # if any remain not found, issue warning
        if(any(is.na(result))) {
            warning("snps not found: ", paste(snp[is.na(result)], collapse=", "))
        }
        return(result)
    }

    # one chromosome:
    snps <- snpinfo$snp
    if(is.null(snps)) snps <- rownames(snpinfo)

    result <- snps[snpinfo$index[match(snp, snps)]]
    if(any(is.na(result))) {
        warning("snps not found: ", paste(snp[is.na(result)], collapse=", "))
    }

    result
}

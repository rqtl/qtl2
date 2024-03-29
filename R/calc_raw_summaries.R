# calculate summary statistics from raw SNP data

#' Calculate estimated heterozygosity from raw SNP genotypes
#'
#' Calculate estimated heterozygosity for each individual from raw SNP genotypes
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param by Indicates whether to summarize by founder strain (`"individual"`) or by marker.
#'
#' @return A vector of heterozygosities, one for each individual or marker.
#'
#' @export
#' @keywords utilities
#' @seealso [recode_snps()], [calc_raw_maf()], [calc_raw_founder_maf()], [calc_raw_geno_freq()], [calc_het()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' DOex_het <- calc_raw_het(DOex)
#' }
calc_raw_het <-
    function(cross, by=c("individual", "marker"))
{
    by <- match.arg(by)

    g <- do.call("cbind", cross$geno)
    g[g < 1 | g > 3] <- NA

    if(by=="individual") return( rowMeans(g==2, na.rm=TRUE) )

    colMeans(g==2, na.rm=TRUE)
}



#' Calculate minor allele frequency from raw SNP genotypes
#'
#' Calculate minor allele frequency from raw SNP genotypes, by individual or by marker
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param by Indicates whether to summarize by founder strain (`"individual"`) or by marker.
#'
#' @return A vector of minor allele frequencies, one for each individual or marker.
#'
#' @export
#' @keywords utilities
#' @seealso [recode_snps()], [calc_raw_founder_maf()], [calc_raw_het()], [calc_raw_geno_freq()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' DOex_maf <- calc_raw_maf(DOex)
#' }
calc_raw_maf <-
    function(cross, by=c("individual", "marker"))
{
    by <- match.arg(by)

    g <- do.call("cbind", cross$geno)
    g[g < 1 | g > 3] <- NA

    if(by=="individual") {
        return( rowMeans(g==3, na.rm=TRUE) + rowMeans(g==2, na.rm=TRUE)/2 )
    }

    # else by marker
    colMeans(g==3, na.rm=TRUE) + colMeans(g==2, na.rm=TRUE)/2

}


#' Calculate founder minor allele frequencies from raw SNP genotypes
#'
#' Calculate minor allele frequency from raw SNP genotypes in founders, by founder strain or by marker
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param by Indicates whether to summarize by founder strain (`"individual"`) or by marker.
#'
#' @return A vector of minor allele frequencies, one for each founder strain or marker.
#'
#' @export
#' @keywords utilities
#' @seealso [recode_snps()], [calc_raw_maf()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' DOex_maf <- calc_raw_founder_maf(DOex)
#' }
calc_raw_founder_maf <-
    function(cross, by=c("individual", "marker"))
{
    if(!("founder_geno" %in% names(cross))) {
        stop("cross does not contain founder genotypes")
    }

    by <- match.arg(by)

    fg <- do.call("cbind", cross$founder_geno)
    fg[fg!=1 & fg!=3] <- NA

    if(by=="individual") {
        return( rowMeans(fg==3, na.rm=TRUE) )
    }

    # else by marker
    colMeans(fg==3, na.rm=TRUE)
}


#' Calculate genotype frequencies from raw SNP genotypes
#'
#' Calculate genotype frequencies from raw SNP genotypes, by individual or by marker
#'
#' @param cross Object of class `"cross2"`. For details, see the
#' [R/qtl2 developer guide](https://kbroman.org/qtl2/assets/vignettes/developer_guide.html).
#' @param by Indicates whether to summarize by individual or by marker.
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#'
#' @return A matrix of genotypes frequencies with 3 columns (AA, AB,
#' and BB) and with rows being either individuals or markers.
#'
#' @export
#' @keywords utilities
#' @seealso [calc_raw_maf()], [calc_raw_het()], [recode_snps()], [calc_geno_freq()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/main/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' gfreq <- calc_raw_geno_freq(DOex)
#' }
calc_raw_geno_freq <-
    function(cross, by=c("individual", "marker"), cores=1)
{
    by <- match.arg(by)

    g <- do.call("cbind", cross$geno)
    g[g<1 | g>3] <- NA

    cores <- setup_cluster(cores)

    if(n_cores(cores)==1) {
        if(by=="individual") {
            result <- t(apply(g, 1, function(a) table(factor(a, levels=1:3))))
        } else {
            result <- t(apply(g, 2, function(a) table(factor(a, levels=1:3))))
        }
    } else {

        if(by=="individual") {

            rn <- rownames(g)
            result <- cluster_lapply(cores, seq_len(nrow(g)), function(row_index) table(factor(g[row_index,], levels=1:3)))
        } else {

            rn <- colnames(g)
            result <- cluster_lapply(cores, seq_len(ncol(g)), function(col_index) table(factor(g[,col_index], levels=1:3)))
        }

        result <- matrix(unlist(result), byrow=TRUE, ncol=3)
        rownames(result) <- rn
    }

    colnames(result) <- c("AA", "AB", "BB")

    result/rowSums(result)
}

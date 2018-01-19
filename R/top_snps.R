#' Create table of top snp associations
#'
#' Create a table of the top snp associations
#'
#' @md
#'
#' @param scan1_output Output of [scan1()].
#' Should contain a component `"snpinfo"`, as when
#' [scan1()] is run with SNP probabilities
#' produced by [genoprob_to_snpprob()].
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
#' @param lodcolumn Selected LOD score column to (a numeric index, or a
#' character string for a column name). Only one value allowed.
#'
#' @param chr Selected chromosome; only one value allowed.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param drop Show all SNPs with LOD score within this amount of the
#' maximum SNP association.
#'
#' @examples
#' \dontrun{
#' # load example DO data from web
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
#' # query function for grabbing info about variants in region
#' dbfile <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
#' query_variants <- create_variant_query_func(dbfile)
#'
#' # SNP association scan, keep information on all SNPs
#' out_snps <- scan1snps(apr, DOex$pmap, DOex$pheno, query_func=query_variants,
#'                       chr=2, start=97, end=98, keep_all_snps=TRUE)
#'
#' # table with top SNPs
#' top_snps(out_snps$lod, out_snps$snpinfo)
#'
#' # top SNPs among the distinct subset at which calculations were performed
#' top_snps(out_snps$lod, out_snps$snpinfo, show_all_snps=FALSE)
#'
#' # top SNPs within 0.5 LOD of max
#' top_snps(out_snps$lod, out_snps$snpinfo, drop=0.5)
#' }
#' @export
#' @seealso [index_snps()], [genoprob_to_snpprob()], [scan1snps()], [plot_snpasso()]
top_snps <-
    function(scan1_output, snpinfo, lodcolumn=1, chr=NULL, drop=1.5, show_all_snps=TRUE)
{
    if(is.null(scan1_output)) stop("scan1_output is NULL")
    if(is.null(snpinfo)) stop("snpinfo is NULL")

    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(length(lodcolumn) > 1) { # If length > 1, take first value
        warning("lodcolumn should have length 1; only first element used.")
        lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(scan1_output))
        if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(scan1_output))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(scan1_output), ")")
    scan1_output <- scan1_output[,lodcolumn,drop=FALSE]

    # reduce to one chromosome
    if(is.null(chr)) {
        uchr <- unique(snpinfo$chr)
        if(length(uchr) > 1) {
            warning("Considering only chr ", uchr[1])
        }
        chr <- uchr[1]
    }
    if(length(chr) > 1) {
        chr <- chr[1]
        warning("Considering only chr ", chr)
    }
    snpinfo <- snpinfo[snpinfo$chr == chr,,drop=FALSE]
    scan1_output <- scan1_output[rownames(scan1_output) %in% snpinfo$snp,,drop=FALSE]

    # check index business
    uindex <- unique(snpinfo$index)
    if(length(uindex) != nrow(scan1_output))
        stop("Something is wrong with snpinfo$index.\n",
             "      length(unique(snpinfo$index)) [",
             length(unique(snpinfo$index)), "] != nrow(scan1_output) [",
             length(scan1_output), "].")

    if(any(snpinfo$index[uindex] != uindex))
        stop("Something is wrong with snpinfo$index.\n",
             "      snpinfo$index[u] should == u for values in snpinfo$index")

    map <- snpinfo_to_map(snpinfo)

    keep <- which(!is.na(scan1_output[,1]) & scan1_output[,1] >= max(scan1_output[ ,1], na.rm=TRUE) - drop)

    # reverse the snp index
    revindex <- rev_snp_index(snpinfo)

   if(show_all_snps) { # expand to all related SNPs
        snpinfo <- snpinfo[revindex %in% keep,,drop=FALSE]
        revindex <- revindex[revindex %in% keep]
        snpinfo$lod <- scan1_output[revindex]
    } else { # just keep the SNPs that were used
        snpinfo$lod <- scan1_output[revindex]
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
    for(i in seq_along(map)) {
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
    index_spl <- split(seq_len(nrow(snpinfo)), snpinfo$index)
    revindex <- rep(seq_along(index_spl), vapply(index_spl, length, 1))
    revindex[unlist(index_spl)] <- revindex

    revindex
}

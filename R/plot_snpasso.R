#' Plot SNP associations
#'
#' Plot SNP associations, with possible expansion from distinct snps to all snps.
#'
#' @param scan1output Output of [scan1()] using
#' SNP probabilities derived by
#' [genoprob_to_snpprob()].
#'
#' @param snpinfo Data frame with SNP information with the following
#'     columns (the last three are generally derived from with
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
#'     groups of SNPs.
#' * `intervals` - Indexes that indicate which marker
#'     intervals the SNPs reside.
#' * `on_map` - Indicate whether SNP coincides with a marker
#'     in the `genoprobs`
#'
#' @param genes Optional data frame containing gene information for
#' the region, with columns `start` and `stop` in Mbp, `strand`
#' (as `"-"`, `"+"`, or `NA`), and `Name`. If included, a
#' two-panel plot is produced, with SNP associations above and
#' gene locations below.
#'
#' @param lodcolumn LOD score column to plot (a numeric index, or a
#' character string for a column name). Only one value allowed.
#'
#' @param show_all_snps If TRUE, expand to show all SNPs.
#'
#' @param chr Vector of character strings with chromosome IDs to plot.
#'
#' @param add If TRUE, add to current plot (must have same map and
#' chromosomes).
#'
#' @param drop_hilit SNPs with LOD score within this amount of the maximum SNP association will be highlighted.
#'
#' @param col_hilit Color of highlighted points
#'
#' @param col Color of other points
#'
#' @param gap Gap between chromosomes. The default is 1% of the total genome length.
#'
#' @param minlod Minimum LOD to display. (Mostly for GWAS, in which
#'     case using `minlod=1` will greatly increase the plotting speed,
#'     since the vast majority of points would be omittted.
#'
#' @param ... Additional graphics parameters.
#'
#' @section Hidden graphics parameters:
#' A number of graphics parameters can be passed via `...`. For
#' example, `bgcolor` to control the background color,`altbgcolor`
#' to control the background color on alternate chromosomes,
#' `altcol` to control the point color on alternate chromosomes,
#' `cex` for character expansion for the points (default 0.5),
#' `pch` for the plotting character for the points (default 16),
#' and `ylim` for y-axis limits.
#'
#' @examples
#' \donttest{
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
#' snp_dbfile <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
#' query_variants <- create_variant_query_func(snp_dbfile)
#'
#' # SNP association scan
#' out_snps <- scan1snps(apr, DOex$pmap, DOex$pheno, query_func=query_variants,
#'                       chr=2, start=97, end=98, keep_all_snps=TRUE)
#'
#' # plot results
#' plot_snpasso(out_snps$lod, out_snps$snpinfo)
#'
#' # can also just type plot()
#' plot(out_snps$lod, out_snps$snpinfo)
#'
#' # plot just subset of distinct SNPs
#' plot(out_snps$lod, out_snps$snpinfo, show_all_snps=FALSE)
#'
#' # highlight the top snps (with LOD within 1.5 of max)
#' plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5)
#'
#' # query function for finding genes in region
#' gene_dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
#' query_genes <- create_gene_query_func(gene_dbfile)
#' genes <- query_genes(2, 97, 98)
#'
#' # plot SNP association results with gene locations
#' plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)
#' }
#'
#' @seealso [plot_scan1()], [plot_coef()], [plot_coefCC()]
#'
#' @importFrom stats setNames
#' @export
#'
plot_snpasso <-
    function(scan1output, snpinfo, genes=NULL, lodcolumn=1, show_all_snps=TRUE, chr=NULL,
             add=FALSE, drop_hilit=NA, col_hilit="violetred", col="darkslateblue",
             gap=NULL, minlod=0, ...)
{
    if(is.null(scan1output)) stop("scan1output is NULL")
    if(is.null(snpinfo)) stop("snpinfo is NULL")
    if(!is_nonneg_number(minlod)) stop("minlod should be a single non-negative number")

    # pull out lod scores
    if(length(lodcolumn)==0) stop("lodcolumn has length 0")
    if(length(lodcolumn) > 1) { # If length > 1, take first value
        warning("lodcolumn should have length 1; only first element used.")
        lodcolumn <- lodcolumn[1]
    }
    if(is.character(lodcolumn)) { # turn column name into integer
        tmp <- match(lodcolumn, colnames(scan1output))
        if(is.na(tmp)) stop('lodcolumn "', lodcolumn, '" not found')
        lodcolumn <- tmp
    }
    if(lodcolumn < 1 || lodcolumn > ncol(scan1output))
        stop("lodcolumn [", lodcolumn, "] out of range (should be in 1, ..., ", ncol(scan1output), ")")
    scan1output <- scan1output[,lodcolumn,drop=FALSE]

    if(!is.null(chr)) { # subset by chr
        if(!any(chr %in% snpinfo$chr)) stop("None of the chr found in snpinfo")
        if(!all(chr %in% snpinfo$chr)) {
            warning("Some chr not found: ", paste(chr[!(chr %in% snpinfo$chr)], collapse=", "))
        }
        snpinfo <- snpinfo[snpinfo$chr %in% chr,,drop=FALSE]
        scan1output <- scan1output[rownames(scan1output) %in% snpinfo$snp,,drop=FALSE]
    }

    if(!is.null(genes)) {
        if(length(unique(snpinfo$chr)) > 1) {
            warning("genes ignored if there are multiple chromosomes")
        } else {
            return( plot_snpasso_and_genes(scan1output, snpinfo, show_all_snps=show_all_snps,
                                           drop_hilit=drop_hilit, col_hilit=col_hilit,
                                           col=col, gap=gap, minlod=minlod, genes=genes, ...) )
        }
    }

    if(nrow(scan1output) == nrow(snpinfo) && all(rownames(scan1output) == snpinfo$snp)) {
        show_all_snps <- FALSE

        snpinfo <- snpinfo[scan1output[,1]>=minlod, , drop=FALSE]
        scan1output <- scan1output[scan1output[,1]>=minlod, 1, drop=FALSE]

        # skip the index business
        # snpinfo -> map
        map <- tapply(seq_len(nrow(snpinfo)), factor(snpinfo$chr, unique(snpinfo$chr)),
                      function(i) setNames(snpinfo$pos[i], snpinfo$snp[i]))
    }
    else {
        snpinfo_spl <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))

        uindex <- unlist(lapply(snpinfo_spl, function(a) unique(a$index)))
        if(length(uindex) != nrow(scan1output)) {
            stop("Something is wrong with snpinfo$index.\n",
                 "      no. unique indexes [",
                 length(uindex), "] != nrow(scan1output) [",
                 nrow(scan1output), "].")
        }

        for(i in seq_along(snpinfo_spl)) {
            uindex <- unique(snpinfo_spl[[i]]$index)
            if(any(snpinfo_spl[[i]]$index[uindex] != uindex)) {
                stop("Something is wrong with snpinfo$index.\n",
                     "      on each chr, index[u] should == u for each index u")
            }
        }

        map <- snpinfo_to_map(snpinfo)

        if(show_all_snps) {
            tmp <- expand_snp_results(scan1output, map, snpinfo)
            scan1output <- tmp$lod
            map <- tmp$map
        }
    }

    # maximum LOD
    maxlod <- max(unclass(scan1output)[,1], na.rm=TRUE)

    # set values < minlod to NA so they're not plotted
    scan1output[scan1output < minlod] <- NA

    # internal function to give defaults to hidden graphics parameters
    plot_snpasso_internal <-
        function(pch=16, cex=0.5, ylim=NULL, bgcolor="gray90",
                 altbgcolor="gray85",
                 drop_hilit=NA, col_hilit="violetred",
                 drop.hilit=NULL, col.hilit=NULL, ...)
    {
        if(!is.null(drop.hilit)) {
            warning("drop.hilit is deprecated; use drop_hilit")
            drop_hilit <- drop.hilit
        }

        if(!is.null(col.hilit)) {
            warning("col.hilit is deprecated; use col_hilit")
            col_hilit <- col.hilit
        }

        if(is.null(ylim))
            ylim <- c(0, maxlod*1.02)

        if(!is.na(drop_hilit) && !is.null(drop_hilit)) {
            col <- c(col, col_hilit)[(scan1output >= maxlod-drop_hilit)+1]
        }

        plot_scan1(scan1output, map, lodcolumn=1, bgcolor=bgcolor, altbgcolor=altbgcolor, ylim=ylim,
                   gap=gap, add=add, col = col, type="p", cex=cex, pch=pch, ...)
    }

    plot_snpasso_internal(drop_hilit=drop_hilit, col_hilit=col_hilit, ...)
}


# expand snp association results according to snpinfo
expand_snp_results <-
    function(snp_results, map, snpinfo)
{
    snpinfo <- split(snpinfo, factor(snpinfo$chr, unique(snpinfo$chr)))

    if(length(map) != length(snpinfo))
        stop("length(map) [", length(map), "] != length(snpinfo) [",
             length(snpinfo), "]")

    if(nrow(snp_results) != length(unlist(map)))
        stop("nrow(snp_results) [", nrow(snp_results), "] != length(unlist(map)) [",
             length(unlist(map)), "]")

    cnames <- rep(names(map), vapply(map, length, 0))
    lodindex <- split(seq_len(nrow(snp_results)), factor(cnames, unique(cnames)))

    result <- NULL
    for(i in seq(along=map)) {
        revindex <- rev_snp_index(snpinfo[[i]])

        map[[i]] <- snpinfo[[i]]$pos
        names(map[[i]]) <- snpinfo[[i]]$snp
        this_result <- unclass(snp_results)[lodindex[[i]],,drop=FALSE][revindex,,drop=FALSE]
        rownames(this_result) <- snpinfo[[i]]$snp

        result <- rbind(result, this_result)
    }

    list(lod=result,
         map=map)
}

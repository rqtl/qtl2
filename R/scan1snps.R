# genome scan at SNPs (or scan of defined region)
#' Single-QTL genome scan at imputed SNPs
#'
#' Perform a single-QTL scan across the genome or a defined region at
#' SNPs genotyped in the founders, by Haley-Knott regression or a
#' liear mixed model, with possible allowance for covariates.
#'
#' @param genoprobs Genotype probabilities as calculated by
#' [calc_genoprob()].
#' @param map Physical map for the positions in the `genoprobs`
#' object: A list of numeric vectors; each vector gives marker
#' positions for a single chromosome.
#' @param pheno A numeric matrix of phenotypes, individuals x phenotypes.
#' @param kinship Optional kinship matrix, or a list of kinship matrices (one
#' per chromosome), in order to use the LOCO (leave one chromosome
#' out) method.
#' @param addcovar An optional numeric matrix of additive covariates.
#' @param Xcovar An optional numeric matrix with additional additive covariates used for
#' null hypothesis when scanning the X chromosome.
#' @param intcovar An optional numeric matrix of interactive covariates.
#' @param weights An optional numeric vector of positive weights for the
#' individuals. As with the other inputs, it must have `names`
#' for individual identifiers.
#' @param reml If `kinship` provided: if `reml=TRUE`, use
#' REML; otherwise maximum likelihood.
#' @param model Indicates whether to use a normal model (least
#' squares) or binary model (logistic regression) for the phenotype.
#' If `model="binary"`, the phenotypes must have values in \eqn{[0, 1]}.
#' @param query_func Function for querying SNP information; see
#' [create_variant_query_func()]). Takes arguments
#' `chr`, `start`, `end`, (with `start` and `end` in the units in
#' `map`, generally Mbp), and returns a data frame containing
#' the columns `snp`, `chr`, `pos`, and `sdp`. (See `snpinfo` below.)
#' @param chr Chromosome or chromosomes to scan
#' @param start Position defining the start of an interval to scan.
#' Should be a single number, and if provided, `chr` should also
#' have length 1.
#' @param end Position defining the end of an interval to scan.
#' Should be a single number, and if provided, `chr` should also
#' have length 1.
#' @param snpinfo Optional data frame of SNPs to scan; if provided,
#' `query_func`, `chr`, `start`, and `end` are ignored. Should
#' contain the following columns:
#' * `chr` - Character string or factor with chromosome
#' * `pos` - Position (in same units as in the `"map"`).
#' * `sdp` - Strain distribution pattern: an integer, between
#'     1 and \eqn{2^n - 2} where \eqn{n} is the number of strains, whose
#'     binary encoding indicates the founder genotypes
#' * `snp` - Character string with SNP identifier (if
#'     missing, the rownames are used).
#' @param batch_length Interval length (in units of `map`, generally
#' Mbp) to scan at one time.
#' @param keep_all_snps SNPs are grouped into equivalence classes based
#'     on position and founder genotypes; if `keep_all_snps=FALSE`,
#'     the return value will contain information only on the indexed
#'     SNPs (one per equivalence class).
#' @param cores Number of CPU cores to use, for parallel calculations.
#' (If `0`, use [parallel::detectCores()].)
#' Alternatively, this can be links to a set of cluster sockets, as
#' produced by [parallel::makeCluster()].
#' @param ... Additional control parameters passed to [scan1()]
#'
#' @return A list with two components: `lod` (matrix of LOD scores)
#' and `snpinfo` (a data frame of SNPs that were scanned,
#' including columns `index` which indicates groups of equivalent
#' SNPs)
#'
#' @details
#' The analysis proceeds as follows:
#' * Call `query_func()` to grab all SNPs over a region.
#' * Use [index_snps()] to group SNPs into equivalence classes.
#' * Use [genoprob_to_snpprob()] to convert `genoprobs` to SNP probabilities.
#' * Use [scan1()] to do a single-QTL scan at the SNPs.
#'
#' @seealso [scan1()], [genoprob_to_snpprob()], [index_snps()], [create_variant_query_func()], [plot_snpasso()]
#'
#' @examples
#' \dontrun{
#' # load example data and calculate genotype probabilities
#' file <- paste0("https://raw.githubusercontent.com/rqtl/",
#'                "qtl2data/master/DOex/DOex.zip")
#' DOex <- read_cross2(file)
#' probs <- calc_genoprob(DOex, error_prob=0.002)
#'
#' snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
#' queryf <- create_variant_query_func(snpdb_file)
#'
#' out <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97, end=98)
#' }
#'
#' @export
#'
scan1snps <-
    function(genoprobs, map, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, model=c("normal", "binary"),
             query_func=NULL, chr=NULL, start=NULL, end=NULL, snpinfo=NULL,
             batch_length=20, keep_all_snps=FALSE, cores=1, ...)
{
    model <- match.arg(model)

    if(is.null(genoprobs)) stop("genoprobs is NULL")
    if(is.null(map)) stop("map is NULL")
    if(is.null(pheno)) stop("pheno is NULL")

    # reduce to common chromosomes
    pchr <- names(genoprobs)
    mchr <- names(map)
    cchr <- pchr[pchr %in% mchr]
    if(length(cchr)==0) {
        warning("No common chromosomes among genoprobs and map")
        return(NULL)
    }
    genoprobs <- genoprobs[,cchr]
    map <- map[cchr]

    # check inputs
    if(length(genoprobs) != length(map)) {
        stop("length(genoprobs) != length(map)")
    }
    if(any(dim(genoprobs)[3,] != vapply(map, length, 1))) {
        stop("genoprobs and map have different numbers of markers")
    }
    dn <- dimnames(genoprobs)[[3]]
    different_names <- FALSE
    for(i in seq_along(dn)) {
        if(any(dn[[i]] != names(map[[i]]))) different_names <- TRUE
    }
    if(different_names) { # different marker names...give a warning (maybe should be an error)
        warning("genoprobs and map have different marker names")
    }
    if(!is.null(snpinfo) && !is.null(query_func)) {
        warning("If snpinfo is provided, chr, start, end, and query_func are all ignored")
    }

    if(!is.null(chr)) chr <- unique(as.character(chr)) # make sure chr is character strings

    if(!is.null(chr) || !is.null(start) || !is.null(end)) { # defined region
        if(!is.null(snpinfo)) warning("If snpinfo provided, chr, start, end, and query_func are all ignored")
        if(is.null(chr)) stop("If start and/or end provided, need to give chr as well")
        if(is.null(query_func)) stop("If chr, start, and/or end provided, need to give query_func too")

        if(length(chr) == 1) {
            if(is.null(start)) start <- -1e15
            if(is.null(end)) end <- 1e15

            if(length(start) > 1 || length(end) > 1) {
                warning("start and end should have length 1; using the first values")
                start <- start[1]
                end <- end[1]
                if(start >= end) stop("Need start < end")
            }

            snpinfo <- query_func(chr, start, end)
        }
    }

    if(length(chr) > 1) {
        if(!is.null(start) || !is.null(end)) {
            warning("If length(chr) > 1, start end end are ignored.")
        }

        if(!all(chr %in% names(map))) {
            stop("Not all chromosomes found: ", paste(chr[!(chr %in% names(map))], collapse=", "))
        }
        map <- map[chr]
        genoprobs <- genoprobs[,chr]
    }

    if(!is.null(snpinfo)) {
        return( scan1snps_snpinfo(genoprobs=genoprobs, map=map, pheno=pheno, kinship=kinship,
                                  addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar,
                                  weights=weights, reml=reml, model=model,
                                  snpinfo=snpinfo, keep_all_snps=keep_all_snps, ...) )
    }

    if(is.null(query_func)) {
        stop("Need to provide either snpinfo or query_func()")
    }

    # split into batches
    if(!is_pos_number(batch_length)) stop("batch_length should be a single positive integer")
    chr <- names(map)
    chr_start <- vapply(map, min, 0)
    chr_end <- vapply(map, max, 0)
    n_batch <- ceiling((chr_end - chr_start)/batch_length)
    endpts <- lapply(seq_along(chr), function(i) seq(chr_start[i], chr_end[i], length.out=n_batch[i]+1))
    batches <- data.frame(chr=rep(chr, n_batch),
                          start=unlist(lapply(endpts, function(a) a[-length(a)])),
                          end=unlist(lapply(endpts, function(a) a[-1])),
                          stringsAsFactors=FALSE)

    # function that does the work
    by_batch_func <- function(i)
    {
        snpinfo <- query_func(batches$chr[i], batches$start[i], batches$end[i])

        if(nrow(snpinfo) == 0) return( list(lod=NULL, snpinfo=snpinfo) )

        # trim off end if necessary
        if(i < nrow(batches) && batches$chr[i] == batches$chr[i+1] &&
           batches$end[i] == batches$start[i+1])
            snpinfo <- snpinfo[snpinfo$pos < batches$end[i], , drop=FALSE]

        scan1snps_snpinfo(genoprobs=genoprobs, map=map, pheno=pheno,
                          kinship=subset_kinship(kinship, chr=batches$chr[i]),
                          addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar,
                          weights=weights, reml=reml, model=model, snpinfo=snpinfo,
                          keep_all_snps=keep_all_snps, ...)
    }

    # set up parallel analysis
    cores <- setup_cluster(cores)

    # the real work
    result <- cluster_lapply(cores, seq_len(nrow(batches)), by_batch_func)

    # combine the multiple results
    omit <- vapply(result, function(a) is.null(a$lod) || is.null(a$snpinfo) || nrow(a$snpinfo)==0, TRUE)
    result <- result[!omit]

    # one useful batch?
    if(length(result)==0) return(NULL)
    if(length(result)==1) return(result[[1]])

    ### lod scores are easy
    lod <- do.call("rbind", lapply(result, "[[", "lod"))

    ### snpinfo: need to revise the index column
    offset <- 0
    for(i in seq_along(result)[-1]) {
        offset <- offset + nrow(result[[i-1]]$snpinfo)
        if(result[[i]]$snpinfo$chr[1] == result[[i-1]]$snpinfo$chr[1]) { # same chromosome so advanced indexes
            result[[i]]$snpinfo$index <- result[[i]]$snpinfo$index + offset
        }
        else offset <- 0
    }
    snpinfo <- do.call("rbind", lapply(result, "[[", "snpinfo"))
    rownames(snpinfo) <- seq_len(nrow(snpinfo))

    list(lod=lod, snpinfo=snpinfo)
}


# scan1snps where snpinfo provided
#     snpinfo must contain snp_id, chr, pos, SDP
scan1snps_snpinfo <-
    function(genoprobs, map, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL, intcovar=NULL,
             weights=NULL, reml=TRUE, model=c("normal", "binary"), snpinfo,
             keep_all_snps=FALSE, ...)
{
    if(is.null(snpinfo) || nrow(snpinfo)==0) return(NULL) # no snps to study
    model <- match.arg(model)

    # reduce to common chromosomes
    pchr <- names(genoprobs)
    mchr <- names(map)
    schr <- unique(snpinfo$chr)
    cchr <- pchr[pchr %in% mchr & pchr %in% schr]
    if(length(cchr)==0) {
        warning("No common chromosomes among genoprobs, map, and snpinfo")
        return(NULL)
    }
    genoprobs <- genoprobs[,cchr]
    map <- map[cchr]
    snpinfo <- snpinfo[snpinfo$chr %in% cchr,,drop=FALSE]

    # snpinfo -> add index
    snpinfo <- index_snps(map, snpinfo)

    # genoprob -> snpprob
    snp_pr <- genoprob_to_snpprob(genoprobs, snpinfo)

    # scan1
    lod <- scan1(snp_pr, pheno=pheno, kinship=subset_kinship(kinship, chr=cchr),
                 addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar,
                 weights=weights, reml=reml, model=model, ...)

    if(!keep_all_snps) {
        snpinfo <- reduce_to_index_snps(snpinfo)
    }

    # return list with lod scores + indexed snpinfo
    list(lod=lod, snpinfo=snpinfo)
}

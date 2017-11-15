# genome scan at SNPs (or scan of defined region)

scan1snps <-
    function(genoprobs, map, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL,
             intcovar=NULL, weights=NULL, reml=TRUE, model=c("normal", "binary"),
             query_func=NULL, chr=NULL, start=NULL, end=NULL, snpinfo=NULL,
             batch_length=20, cores=1, keep_all_snps=FALSE, ...)
{
    model <- match.arg(model)

    # check inputs
    if(!is.null(snpinfo) && !is.null(query_func)) {
        warning("If snpinfo is provided, chr, start, end, and query_func are all ignored")
    }

    if(!is.null(chr) || !is.null(start) || !is.null(end)) { # defined region
        if(!is.null(snpinfo)) warning("If snpinfo provided, chr, start, end, and query_func are all ignored")
        if(is.null(chr)) stop("If start and/or end provided, need to give chr as well")
        if(is.null(query_func)) stop("If chr, start, and/or end provided, need to give query_func too")
        if(is.null(start)) start <- -1e15
        if(is.null(end)) end <- 1e15

        snpinfo <- query_func(chr, start, end)
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
    if(batch_length <= 0) stop("batch_length should be > 0")
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

        scan1snps_snpinfo(genoprobs=genoprobs, map=map, pheno=pheno, kinship=kinship,
                          addcovar=addcovar, Xcovar=Xcovar, intcovar=intcovar,
                          weights=weights, reml=reml, model=model, snpinfo=snpinfo, ...)
    }

    # set up parallel analysis
    cores <- setup_cluster(cores)

    result <- cluster_lapply(cores, 1:nrow(batches), by_batch_func)

    # was one batch?
    if(length(result)==1) return(result)

    # combine the multiple results
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
    rownames(snpinfo) <- 1:nrow(snpinfo)

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

    # snpinfo -> add index
    snpinfo <- index_snps(map, snpinfo)

    # genoprob -> snpprob
    snp_pr <- genoprob_to_snpprob(genoprobs, snpinfo)

    # scan1
    lod <- scan1(snp_pr, pheno=pheno, kinship=kinship, addcovar=addcovar, Xcovar=Xcovar,
                 intcovar=intcovar, weights=weights, reml=reml, model=model, ...)

    if(!keep_all_snps) {
        snpinfo <- snpinfo[unique(snpinfo$index),,drop=FALSE]
        snpinfo$index <- 1:nrow(snpinfo)
    }

    # return list with lod scores + indexed snpinfo
    list(lod=lod, snpinfo=snpinfo)
}
